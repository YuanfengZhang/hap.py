#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Derive a versioned summary for the fixed GRCh38 fixture using validated insertion-boundary rules."""
import argparse
import csv
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path
import pysam


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def project(path: Path):
    """Yield explicit record projections, rejecting duplicate identities."""
    seen: set[tuple] = set()
    with pysam.VariantFile(str(path)) as source:
        for record in source:
            key: tuple = (record.contig, record.pos, record.ref, tuple(record.alts or ()))
            if key in seen:
                raise AssertionError(f'Duplicate record: {key}')
            seen.add(key)
            samples: dict = {sample: {field: record.samples[sample].get(field)
                            for field in ('GT', 'BD', 'BK', 'BVT', 'BLT')}
                            for sample in ('TRUTH', 'QUERY')}
            yield json.loads(json.dumps([key, tuple(record.filter), samples]))


def compare_records(actual: Path, expected: Path) -> int:
    """Compare ordered fields against the independently verified legacy baseline."""
    from itertools import zip_longest

    count: int = 0
    with gzip.open(expected, 'rt') as stream:
        for index, (row, line) in enumerate(zip_longest(project(actual), stream), 1):
            if row is None or line is None:
                raise AssertionError(f'Unexpected record or classification change: record count at {index}')
            wanted: list = json.loads(line)
            for field, value, reference in zip(('variant key', 'FILTER', 'samples'), row, wanted):
                if value != reference:
                    raise AssertionError(f'Unexpected record or classification change: record {index} '
                                         f'{field}: expected {reference!r}, observed {value!r}')
            count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--actual-vcf', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root: Path = Path(__file__).resolve().parents[2]
    fixture: Path = root / 'example/happy'
    contract: dict = json.loads((root / 'src/data/insertion-boundary-summary-v1.json').read_text())
    for name, checksum in contract['input_sha256'].items():
        if not (digest(fixture / name) == checksum):
            raise AssertionError('Fixed input changed: ' + name)
    intervals: dict[str, list[tuple[int, int]]] = {}
    with gzip.open(fixture / 'PG_Conf_hg38-chr21.bed.gz', 'rt') as stream:
        for line in stream:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            chrom, start, end, *_ = line.split()
            intervals.setdefault(chrom, []).append((int(start), int(end)))
    wanted: set[tuple] = set()
    for chrom, position, ref, alts in contract['insertion_keys']:
        key: tuple = (chrom, position, ref, tuple(alts))
        if not (len(ref) == 1 and all(alt.startswith(ref) and len(alt) > 1 for alt in alts)):
            raise AssertionError(key)
        start, end = position - 1, position + 1
        covered: list[bool] = [any(a <= base < b for a, b in intervals[chrom])
                               for base in range(start, end)]
        if not (any(covered) and not all(covered)):
            raise AssertionError((key, covered))
        wanted.add(key)
    count = compare_records(args.actual_vcf, root / 'src/data' / contract['expected_records'])
    if not (count == contract['record_count']):
        raise AssertionError('Unexpected record count')
    deltas: dict[str, int] = {'ALL': 0, 'PASS': 0}
    with pysam.VariantFile(str(args.actual_vcf)) as source:
        found: set[tuple] = set()
        for record in source:
            key = (record.contig, record.pos, record.ref, tuple(record.alts or ()))
            if key not in wanted:
                continue
            query = record.samples['QUERY']
            if not (query['BD'] == 'UNK' and query['BK'] == 'lm' and query['BVT'] == 'INDEL'):
                raise AssertionError(key)
            if not (query['GT'] and any(allele and allele > 0 for allele in query['GT'])):
                raise AssertionError(key)
            if not ('CONF' not in record.info.get('Regions', ())):
                raise AssertionError(key)
            found.add(key)
            deltas['ALL'] += 1
            if not record.filter or tuple(record.filter) == ('PASS',):
                deltas['PASS'] += 1
        if not (found == wanted):
            raise AssertionError('Missing boundary record')
    with (fixture / 'expected-stratified.summary.csv').open() as stream:
        reader = csv.DictReader(stream)
        fields: list[str] = list(reader.fieldnames or [])
        rows: list[dict] = list(reader)
    for row in rows:
        if row['Type'] != 'INDEL':
            continue
        delta: int = deltas[row['Filter']]
        old_fp, old_unk, total = (Decimal(row[key]) for key in ('QUERY.FP', 'QUERY.UNK', 'QUERY.TOTAL'))
        query_tp: Decimal = total - old_fp - old_unk
        fp, unk = old_fp - delta, old_unk + delta
        if not (fp >= 0):
            raise AssertionError('Negative false-positive count')
        row['QUERY.FP'], row['QUERY.UNK'] = str(fp), str(unk)
        precision: Decimal = query_tp / (query_tp + fp) if query_tp + fp else Decimal(0)
        truth_tp, truth_fn = Decimal(row['TRUTH.TP']), Decimal(row['TRUTH.FN'])
        recall: Decimal = truth_tp / (truth_tp + truth_fn) if truth_tp + truth_fn else Decimal(0)
        row['METRIC.Precision'] = str(precision)
        row['METRIC.Frac_NA'] = str(unk / total if total else Decimal(0))
        row['METRIC.F1_Score'] = str(2 * precision * recall / (precision + recall) if precision + recall else Decimal(0))
    with args.output.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    receipt: dict = {'contract': contract['version'], 'comparison': 'ordered explicit record fields',
                     'record_count': count, 'derived_reclassifications': deltas,
                     'golden_unchanged': True, 'scope': 'Fixed-input summary only; extended metrics remain independently and strictly validated'}
    Path(str(args.output) + '.receipt.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')


if __name__ == '__main__':
    main()
