#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate gene-aggregate and exon labels for numeric BED suffixes."""
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
    vcf: Path = work / 'input.vcf'
    header: str = '##fileformat=VCFv4.2\n##contig=<ID=chrT,length=100>\n'
    for field in ('GT', 'BD', 'BK'):
        header += f'##FORMAT=<ID={field},Number=1,Type=String,Description="{field}">\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n'
    vcf.write_text(header + 'chrT\t10\t.\tC\tT\t40\tPASS\t.\tGT:BD:BK\t0/1:TP:gm\t0/1:TP:gm\n')
    compressed: Path = work / 'input.vcf.gz'
    bed: Path = work / 'regions.bed'
    labels: list[str] = ['GENE_1', 'GENE_12', 'PLAIN', 'GENE2', '123', '_1', 'GENE_']
    bed.write_text(''.join(f'chrT\t0\t20\t{label}\n' for label in labels))
    table: Path = work / 'counts.tsv'
    commands: list[tuple[str, list[str]]] = [
        ('compress', [str(binary / 'bcftools'), 'view', '-Oz', '-o', str(compressed), str(vcf)]),
        ('index', [str(binary / 'bcftools'), 'index', '-t', str(compressed)]),
        ('quantify', [str(binary / 'quantify'), str(compressed), '-r', str(reference),
                      '--type', 'ga4gh', '-R', 'exons:' + str(bed), '-o', str(table)]),
    ]
    results: list[dict] = []
    for step, command in commands:
        with (work / (step + '.log')).open('w') as log:
            code: int = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False).returncode
        results.append({'step': step, 'command': command, 'exit_code': code})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        if code:
            raise SystemExit(code)
    with table.open() as stream:
        rows: dict[str, dict] = {row['Subset']: row for row in csv.DictReader(stream, delimiter='\t')
                                 if row['Type'] == 'SNP' and row['Filter'] == 'ALL'
                                 and all(row[key] == '*' for key in ('Subtype', 'Genotype', 'QQ'))}
    expected: set[str] = {'*', 'exons', 'exons_GENE'} | {'exons_' + label for label in labels}
    if not (set(rows) == expected):
        raise AssertionError((set(rows), expected))
    for label, row in rows.items():
        if not (Decimal(row['TRUTH.TP']) == 1 and Decimal(row['QUERY.TP']) == 1):
            raise AssertionError((label, row))
        if label.startswith('exons_'):
            expected_level: int = 2 if label in ('exons_GENE_1', 'exons_GENE_12') else 1
            if not (Decimal(row['Subset.Level']) == expected_level):
                raise AssertionError((label, row['Subset.Level']))
    # Region size sums input interval lengths; do not confuse this with deduplicated variant counts.
    for label, row in rows.items():
        if label == '*':
            continue
        expected_size: int = 140 if label == 'exons' else 40 if label == 'exons_GENE' else 20
        if not (Decimal(row['Subset.Size']) == expected_size):
            raise AssertionError((label, row['Subset.Size']))
    results.append({'label_contract_pass': True, 'labels': sorted(rows)})
    (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    print('BED region hierarchy: aggregate, exon and unsplit controls PASS')


if __name__ == '__main__':
    main()
