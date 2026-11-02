#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate quantified FP.al and FP.gt definitions using known BD/BK classifications."""
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
    for field, label in (('GT', 'Genotype'), ('BD', 'Decision'), ('BK', 'Kind')):
        header += f'##FORMAT=<ID={field},Number=1,Type=String,Description="{label}">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n'
    results: list[dict] = []
    cases: list[tuple[str, str, str, tuple[int, int, int]]] = [
        ('allele', 'FP', 'lm', (1, 0, 1)),
        ('genotype', 'FP', 'am', (0, 1, 1)),
        ('match', 'TP', 'lm', (0, 0, 0)),
        ('plain', 'FP', '.', (0, 0, 1)),
    ]
    for name, decision, kind, expected in cases:
        vcf: Path = work / (name + '.vcf')
        compressed: Path = work / (name + '.vcf.gz')
        table: Path = work / (name + '.tsv')
        vcf.write_text(header + 'chrT\t2\t.\tC\tT\t40\tPASS\t.\tGT:BD:BK\t'
                       + f'0/1:TP:gm\t0/1:{decision}:{kind}\n')
        commands: list[tuple[str, list[str]]] = [
            ('compress', [str(binary / 'bcftools'), 'view', '-Oz', '-o', str(compressed), str(vcf)]),
            ('index', [str(binary / 'bcftools'), 'index', '-t', str(compressed)]),
            ('quantify', [str(binary / 'quantify'), str(compressed), '-r', str(reference),
                          '--type', 'ga4gh', '-o', str(table)]),
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
                                if row['Type'] == 'SNP' and row['Filter'] == 'ALL'
                                and all(row[key] == '*' for key in ('Subtype', 'Subset', 'Genotype', 'QQ'))]
        if not (len(rows) == 1):
            raise AssertionError((name, len(rows)))
        actual: tuple[Decimal, ...] = tuple(Decimal(rows[0][key])
                                            for key in ('FP.al', 'FP.gt', 'QUERY.FP'))
        if not (actual == expected):
            raise AssertionError((name, actual, expected))
        results.append({'case': name, 'expected': expected, 'classification_pass': True})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    print('FP classification: four controls PASS')


if __name__ == '__main__':
    main()
