#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate three installed CLI engines, thread consistency, and independently calculated metrics."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess

import pysam


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, default=Path('/opt/happy/bin'))
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    work: Path = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    reference: Path = work / 'reference.fa'
    reference.write_text('>custom\n' + 'ACGT' * 100 + '\n')
    pysam.faidx(str(reference))
    header: str = ('##fileformat=VCFv4.2\n##contig=<ID=custom,length=400>\n'
                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
    truth: Path = work / 'truth.vcf'
    query: Path = work / 'query.vcf'
    truth.write_text(header + 'custom\t2\t.\tC\tT\t40\tPASS\t.\tGT\t0/1\n'
                     + 'custom\t10\t.\tC\tA\t20\tPASS\t.\tGT\t0/1\n')
    query.write_text(header + 'custom\t2\t.\tC\tT\t40\tPASS\t.\tGT\t0/1\n'
                     + 'custom\t18\t.\tC\tA\t10\tPASS\t.\tGT\t0/1\n')
    bed: Path = work / 'confident.bed'
    bed.write_text('custom\t0\t400\n')
    results: list[dict] = []
    for engine in ('xcmp', 'vcfeval', 'scmp-somatic'):
        expected_table: list[dict] | None = None
        for threads in (1, 2):
            prefix: Path = work / f'{engine}-{threads}'
            command: list[str] = [str(args.bin_dir.resolve() / 'hap.py'), str(truth), str(query),
                '-r', str(reference), '-f', str(bed), '-o', str(prefix), '--engine', engine,
                '--threads', str(threads), '--scratch-prefix', str(work)]
            with prefix.with_suffix('.log').open('w') as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            with Path(str(prefix) + '.summary.csv').open() as stream:
                rows: list[dict] = list(csv.DictReader(stream))
            row: dict = next(item for item in rows if item['Type'] == 'SNP' and item['Filter'] == 'ALL')
            expected: dict[str, float] = {'TRUTH.TOTAL': 2, 'TRUTH.TP': 1, 'TRUTH.FN': 1,
                'QUERY.TOTAL': 2, 'QUERY.FP': 1, 'QUERY.UNK': 0,
                'METRIC.Recall': 0.5, 'METRIC.Precision': 0.5, 'METRIC.F1_Score': 0.5}
            for field, value in expected.items():
                if float(row[field]) != value:
                    raise ValueError(f'{engine}/{threads}: {field}={row[field]}, expected {value}')
            for suffix in ('.extended.csv', '.roc.all.csv.gz', '.vcf.gz'):
                path: Path = Path(str(prefix) + suffix)
                if not path.is_file() or path.stat().st_size == 0:
                    raise ValueError(f'Missing output: {path}')
            with Path(str(prefix) + '.extended.csv').open() as stream:
                extended: list[dict] = list(csv.DictReader(stream))
            aggregate: dict = next(item for item in extended if item['Type'] == 'SNP'
                and item['Subtype'] == '*' and item['Genotype'] == '*'
                and item['Subset'] == '*' and item['Filter'] == 'ALL')
            if float(aggregate['QUERY.TP']) != 1:
                raise ValueError(f'{engine}/{threads}: wrong extended QUERY.TP')
            if expected_table is not None and rows != expected_table:
                raise ValueError(f'{engine}: thread-dependent summary')
            expected_table = rows
            results.append({'engine': engine, 'threads': threads, 'exit_code': 0})
            (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
            print(f'PASS {engine} threads={threads}', flush=True)


if __name__ == '__main__':
    main()
