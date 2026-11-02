#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate known answers requiring CONF coverage of both reference flanks of pure insertions."""
import argparse
import csv
from decimal import Decimal
import json
from pathlib import Path
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    binary: Path = args.bin_dir.resolve()
    work: Path = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    reference: Path = work / 'reference.fa'
    reference.write_text('>chrT\n' + 'ACGT' * 25 + '\n')
    header: str = '##fileformat=VCFv4.2\n##contig=<ID=chrT,length=100>\n'
    for field in ('GT', 'BD', 'BK'):
        header += f'##FORMAT=<ID={field},Number=1,Type=String,Description="{field}">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n'
    results: list[dict] = []
    cases: list[tuple[str, str, int, int, str, tuple[int, int]]] = [
        ('left', 'CT', 10, 20, 'INDEL', (0, 1)),
        ('right', 'CT', 0, 10, 'INDEL', (0, 1)),
        ('contained', 'CT', 9, 11, 'INDEL', (1, 0)),
        ('snp', 'T', 9, 10, 'SNP', (1, 0)),
    ]
    for name, alt, start, end, variant_type, expected in cases:
        vcf: Path = work / (name + '.vcf')
        compressed: Path = work / (name + '.vcf.gz')
        table: Path = work / (name + '.tsv')
        bed: Path = work / (name + '.bed')
        bed.write_text(f'chrT\t{start}\t{end}\n')
        vcf.write_text(header + f'chrT\t10\t.\tC\t{alt}\t40\tPASS\t.\tGT:BD:BK\t./.:.:.\t0/1:FP:lm\n')
        commands: list[tuple[str, list[str]]] = [
            ('compress', [str(binary / 'bcftools'), 'view', '-Oz', '-o', str(compressed), str(vcf)]),
            ('index', [str(binary / 'bcftools'), 'index', '-t', str(compressed)]),
            ('quantify', [str(binary / 'quantify'), str(compressed), '-r', str(reference),
                          '--type', 'ga4gh', '-R', 'CONF:' + str(bed), '-o', str(table)]),
        ]
        for step, command in commands:
            with (work / f'{name}-{step}.log').open('w') as log:
                code: int = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                           check=False).returncode
            results.append({'case': name, 'step': step, 'command': command, 'exit_code': code})
            (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
            if code:
                raise SystemExit(code)
        with table.open() as stream:
            rows: list[dict] = [row for row in csv.DictReader(stream, delimiter='\t')
                                if row['Type'] == variant_type and row['Filter'] == 'ALL'
                                and all(row[key] == '*' for key in ('Subtype', 'Subset', 'Genotype', 'QQ'))]
        if not (len(rows) == 1):
            raise AssertionError((name, len(rows)))
        actual: tuple[Decimal, ...] = tuple(Decimal(rows[0][key]) for key in ('QUERY.FP', 'QUERY.UNK'))
        if not (actual == expected):
            raise AssertionError((name, actual, expected))
        results.append({'case': name, 'expected_fp_unk': expected, 'boundary_pass': True})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    print('Insertion CONF boundary: four controls PASS')


if __name__ == '__main__':
    main()
