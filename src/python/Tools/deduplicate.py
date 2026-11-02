# -*- coding: utf-8 -*-
"""Remove complete duplicate records with memory bounded by one coordinate."""
from contextlib import ExitStack
import gzip
import logging
import os
from pathlib import Path
import tempfile
from typing import Iterable, Iterator

import pysam

from Tools.bcftools import runBcftools
from Tools.scratch import ScratchWorkspace


def unique_records(lines: Iterable[str]) -> Iterator[str]:
    """Keep first full record at each coordinate; reject coordinate regressions."""
    previous: tuple[str, int] | None = None
    finished: set[str] = set()
    seen: set[str] = set()
    for line in lines:
        fields = line.split('\t', 2)
        coordinate = (fields[0], int(fields[1]))
        if coordinate != previous:
            if previous is not None:
                if coordinate[0] == previous[0]:
                    if coordinate[1] < previous[1]:
                        raise ValueError('Duplicate removal requires coordinate-sorted records')
                else:
                    finished.add(previous[0])
                    if coordinate[0] in finished:
                        raise ValueError('Duplicate removal requires contiguous chromosomes')
            previous = coordinate
            seen.clear()
        key = line.rstrip('\r\n')
        if key not in seen:
            seen.add(key)
            yield line


def deduplicate(filename: str) -> tuple[int, int]:
    """Replace output only after complete writing/indexing; return total/removed.

    VCF text/header is preserved. BCF uses HTSlib's VCF serialization as its key.
    No global sorting or allele-only equality is performed. Data and index
    replacement are separate operations, not a multi-file atomic transaction.
    """
    path = Path(filename)
    is_bcf = filename.endswith('.bcf')
    compressed = filename.endswith('.gz')
    suffix = '.bcf' if is_bcf else '.vcf.gz' if compressed else '.vcf'
    fd, temporary = tempfile.mkstemp(prefix='.deduplicate-', suffix=suffix, dir=path.parent)
    os.close(fd)
    total = retained = 0
    with ScratchWorkspace() as scratch:
        scratch.register(temporary)
        for extension in ('.csi', '.tbi'):
            scratch.register(temporary + extension)
        with ExitStack() as stack:
            if is_bcf:
                source = stack.enter_context(pysam.VariantFile(filename))
                dest = stack.enter_context(pysam.VariantFile(temporary, 'wb', header=source.header))
                # Keep the native object for writing while comparing its full text.
                previous = None
                finished: set[str] = set()
                seen: set[str] = set()
                for record in source:
                    total += 1
                    coordinate = (record.contig, record.pos)
                    if coordinate != previous:
                        if previous:
                            if coordinate[0] == previous[0] and coordinate[1] < previous[1]:
                                raise ValueError('Duplicate removal requires coordinate-sorted records')
                            if coordinate[0] != previous[0]:
                                finished.add(previous[0])
                                if coordinate[0] in finished:
                                    raise ValueError('Duplicate removal requires contiguous chromosomes')
                        seen.clear()
                        previous = coordinate
                    key = str(record).rstrip('\r\n')
                    if key not in seen:
                        seen.add(key)
                        dest.write(record)
                        retained += 1
            else:
                source = stack.enter_context(gzip.open(filename, 'rt', encoding='utf-8', newline='')
                                             if compressed else open(filename, encoding='utf-8', newline=''))
                dest = stack.enter_context(pysam.BGZFile(temporary, 'w') if compressed
                                           else open(temporary, 'wb'))

                def records() -> Iterator[str]:
                    nonlocal total
                    for line in source:
                        if line.startswith('#'):
                            dest.write(line.encode())
                        else:
                            total += 1
                            yield line
                for line in unique_records(records()):
                    dest.write(line.encode())
                    retained += 1
        if total != retained:
            index_suffix = '.csi' if is_bcf else '.tbi' if compressed else None
            if index_suffix:
                runBcftools('index', *([] if is_bcf else ['-t']), temporary)
            os.replace(temporary, filename)
            for extension in ('.tbi', '.csi'):
                Path(filename + extension).unlink(missing_ok=True)
            if index_suffix:
                os.replace(temporary + index_suffix, filename + index_suffix)
    logging.info('Complete-record deduplication %s: total=%d retained=%d removed=%d',
                 filename, total, retained, total - retained)
    return total, total - retained
