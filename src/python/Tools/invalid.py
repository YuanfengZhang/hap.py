# -*- coding: utf-8 -*-
"""Stream invalid GT/AD records into durable, auditable quarantine files."""
from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass
import gzip
import logging
import os
from pathlib import Path
import re
import shutil
import tempfile
from typing import Iterator

import pysam


@dataclass(frozen=True)
class SplitResult:
    valid_path: str
    total: int
    invalid: int


def check_outputs(outputs: list[str], protected: list[str]) -> None:
    """Reject collisions and existing artifacts before processing any input."""
    raw_paths = [Path(path) for output in outputs
                 for path in (output, output + '.reasons.tsv')]
    paths = [path.resolve() for path in raw_paths]
    if len(set(paths)) != len(paths) or set(paths) & {Path(p).resolve() for p in protected}:
        raise ValueError('Invalid-record output collides with an input or output path')
    for path in raw_paths:
        if path.exists() or path.is_symlink():
            raise FileExistsError(f'Refusing to overwrite invalid-record output: {path}')
        if not path.parent.is_dir():
            raise FileNotFoundError(f'Invalid-record output directory does not exist: {path.parent}')


def record_errors(line: str, samples: list[str]) -> list[tuple[str, str, str]]:
    """Validate raw fields, never use pysam GT tuples (invalid indexes become None)."""
    fields = line.rstrip('\r\n').split('\t')
    expected = 9 + len(samples) if samples else 8
    if len(fields) != expected or int(fields[1]) < 1:
        raise ValueError('Malformed VCF record: field count or position')
    alleles = 1 if fields[4] == '.' else 1 + len(fields[4].split(','))
    if not samples:
        return []
    formats = fields[8].split(':')
    errors: list[tuple[str, str, str]] = []
    for name, value in zip(samples, fields[9:], strict=True):
        values = value.split(':')
        if len(values) > len(formats):
            raise ValueError('Malformed VCF sample: more values than FORMAT fields')
        data = dict(zip(formats, values))
        gt = data.get('GT', '.')
        genotype = re.split(r'[/|]', gt)
        if any(g != '.' and re.fullmatch(r'[0-9]+', g) is None for g in genotype):
            raise ValueError(f'Malformed GT syntax: {gt}')
        if len(genotype) > 2:  # BCFHelpers.hh MAX_GT
            errors.append((name, 'GT_PLOIDY', f'ploidy={len(genotype)} capacity=2'))
        if any(g != '.' and int(g) >= alleles for g in genotype):
            errors.append((name, 'GT_ALLELE', f'GT={gt} allele_count={alleles}'))
        ad = data.get('AD', '.')
        if ad != '.':
            depths = ad.split(',')
            if any(d != '.' and re.fullmatch(r'-?[0-9]+', d) is None for d in depths):
                raise ValueError(f'Malformed AD syntax: {ad}')
            if len(depths) > alleles:
                errors.append((name, 'AD_CAPACITY', f'values={len(depths)} allele_count={alleles}'))
    return errors


def split_invalid(source: str, output: str, valid_path: str) -> SplitResult:
    """Split readable VCF/BCF; only publish quarantine after complete parsing.

    VCF quarantine preserves raw header/record text. BCF is serialized through
    HTSlib. Parsing/I/O failures remain fatal; no partially scanned final output.
    The caller owns valid_path and must clean it, including on failure.
    """
    check_outputs([output], [source, valid_path])
    pending: list[str] = []
    counts: Counter[str] = Counter()
    total = invalid = 0
    try:
        with ExitStack() as stack:
            with open(source, 'rb') as probe:
                magic = probe.read(18)
            compressed = magic[:2] == b'\x1f\x8b'
            parse_source = source
            if compressed and magic[12:14] != b'BC':
                fd, parse_source = tempfile.mkstemp(prefix='invalid-reader-', suffix='.vcf.gz',
                                                   dir=str(Path(valid_path).parent))
                os.close(fd)
                pending.append(parse_source)
                with gzip.open(source, 'rb') as src, pysam.BGZFile(parse_source, 'w') as dest:
                    shutil.copyfileobj(src, dest)
            variants = stack.enter_context(pysam.VariantFile(parse_source))
            samples = list(variants.header.samples)
            contigs = set(variants.header.contigs)
            declared_formats = set(variants.header.formats)
            declared_infos = set(variants.header.info)
            declared_filters = set(variants.header.filters)

            def checked_records() -> Iterator[pysam.VariantRecord]:
                for record in variants:
                    if record.contig not in contigs:
                        raise ValueError('Record contig is absent from header')
                    if (set(record.format) - declared_formats or
                            set(record.info) - declared_infos or
                            set(record.filter) - declared_filters):
                        raise ValueError('Record FORMAT/INFO/FILTER is absent from header')
                    yield record
            if variants.is_bcf:
                header = str(variants.header)
                lines: Iterator[str] = (str(record) for record in checked_records())
            else:
                raw = stack.enter_context(gzip.open(source, 'rt', encoding='utf-8', newline='')
                                          if compressed else open(source, encoding='utf-8', newline=''))
                headers: list[str] = []
                for line in raw:
                    if not line.startswith('#'):
                        raise ValueError('VCF missing #CHROM header')
                    headers.append(line)
                    if line.startswith('#CHROM\t'):
                        break
                else:
                    raise ValueError('VCF missing #CHROM header')
                header = ''.join(headers)

                def checked_lines() -> Iterator[str]:
                    for line, record in zip(raw, checked_records(), strict=True):
                        # Force HTSlib parsing without its lossy GT tuple conversion.
                        yield line
                lines = checked_lines()
            artifacts: list[str] = []
            for suffix in ('.vcf.gz', '.tsv'):
                fd, name = tempfile.mkstemp(prefix='.invalid-pending-', suffix=suffix,
                                            dir=str(Path(output).parent))
                os.close(fd)
                pending.append(name)
                artifacts.append(name)
            rejected = stack.enter_context(pysam.BGZFile(artifacts[0], 'w'))
            reasons = stack.enter_context(open(artifacts[1], 'w', encoding='utf-8', newline=''))
            valid = stack.enter_context(pysam.BGZFile(valid_path, 'w'))
            valid.write(header.encode('utf-8'))
            rejected.write(header.encode('utf-8'))
            reasons.write('record\tCHROM\tPOS\tsample\treason\tdetail\n')
            for total, line in enumerate(lines, 1):
                errors = record_errors(line, samples)
                if errors:
                    invalid += 1
                    rejected.write(line.encode('utf-8'))
                    fields = line.split('\t', 2)
                    for sample, reason, detail in errors:
                        counts[reason] += 1
                        reasons.write(f'{total}\t{fields[0]}\t{fields[1]}\t{sample}\t{reason}\t{detail}\n')
                else:
                    valid.write(line.encode('utf-8'))
        # Hard links publish without replacing files created by a concurrent run.
        os.link(artifacts[0], output)
        try:
            os.link(artifacts[1], output + '.reasons.tsv')
        except OSError:
            os.unlink(output)
            raise
        logging.warning('Invalid-record scan %s: total=%d valid=%d excluded=%d reasons=%s; '
                        'evaluation excludes these records; quarantine=%s',
                        source, total, total - invalid, invalid, dict(counts), output)
        if total == invalid:
            raise ValueError(f'No valid records remain after invalid-record scan: {source}; '
                             f'quarantine preserved at {output}')
        return SplitResult(valid_path, total, invalid)
    finally:
        for path in pending:
            Path(path).unlink(missing_ok=True)
