#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate CONF union v1: coverage on both insertion flanks, two-sample classification, and BS propagation."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from decimal import Decimal
import subprocess
import pysam


def check_records(path: Path, fixed: bool) -> None:
    expected_identity = [(11, 'C', (0, 1), (0, 1)), (19, 'C', (0, 1), (0, 1)),
                         (20, 'ATT', (0, 1), (0, 1)), (21, 'ATT', (0, 1), (None, None)),
                         (22, 'C', (None, None), (0, 1))]
    with pysam.VariantFile(str(path)) as stream:
        records = list(stream)
    identity = [(r.pos, r.alts[0], r.samples['TRUTH'].get('GT'), r.samples['QUERY'].get('GT')) for r in records]
    if identity != expected_identity or any(r.contig != 'chrT' or r.ref != 'A' or len(r.alts) != 1 for r in records):
        raise ValueError('Record identity or count changed')
    for record in records:
        confident = fixed or record.pos <= 20
        expected_regions = {'TS_contained' if fixed else 'TS_boundary'}
        if confident:
            expected_regions.add('CONF')
        if set(record.info.get('Regions', ())) != expected_regions:
            raise ValueError(f'{record.pos}: Regions differ from the CONF-union-v1 contract')
        expected = ('TP', 'TP') if record.pos <= 20 else (
            ('FN' if fixed else 'UNK', '.') if record.pos == 21 else
            ('.', 'FP' if fixed else 'UNK'))
        actual = tuple(record.samples[name].get('BD') for name in ('TRUTH', 'QUERY'))
        if actual != expected:
            raise ValueError(f'{record.pos}: BD {actual!r} != {expected!r}')


def check_counts(rows: list[dict[str, str]], fixed: bool) -> None:
    candidates = [row for row in rows if all(row[k] == '*' for k in ('Subset', 'Subtype', 'Genotype', 'QQ')) and row['Filter'] == 'ALL']
    selected = {row['Type']: row for row in candidates}
    if len(candidates) != len(selected):
        raise ValueError('Duplicate count type')
    expected = {
        'SNP': {'TRUTH.TOTAL': 2, 'TRUTH.TP': 2, 'TRUTH.FN': 0,
                'QUERY.TOTAL': 3, 'QUERY.TP': 2, 'QUERY.FP': int(fixed), 'QUERY.UNK': int(not fixed), 'FP.al': int(fixed)},
        'INDEL': {'TRUTH.TOTAL': 2 if fixed else 1, 'TRUTH.TP': 1, 'TRUTH.FN': int(fixed),
                  'QUERY.TOTAL': 1, 'QUERY.TP': 1, 'QUERY.FP': 0, 'QUERY.UNK': 0}}
    if selected.keys() != expected.keys():
        raise ValueError('Missing or duplicate count type')
    for typ, columns in expected.items():
        for column, value in columns.items():
            if Decimal(selected[typ][column]) != value:
                raise ValueError(f'{typ}/{column}: Counts violate CONF-union-v1')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', required=True, type=Path)
    parser.add_argument('--work-dir', required=True, type=Path)
    args = parser.parse_args()
    tools = args.bin_dir.resolve()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    results: list[dict] = []

    def run(name: str, command: list[str]) -> None:
        process = subprocess.run(command, cwd=work, capture_output=True, text=True, check=False)
        (work / f'{name}.stdout').write_text(process.stdout)
        (work / f'{name}.stderr').write_text(process.stderr)
        results.append({'name': name, 'command': command, 'exit_code': process.returncode})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        if process.returncode:
            raise RuntimeError(f'{name} exited {process.returncode}: {process.stderr}')

    (work / 'ref.fa').write_text('>chrT\n' + 'A' * 128 + '\n')
    pysam.faidx(str(work / 'ref.fa'))
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chrT,length=128>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    (work / 'truth.vcf').write_text(
        header + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n'
        + 'chrT\t20\t.\tA\tATT\t30\tPASS\t.\tGT\t0/1\n'
        + 'chrT\t21\t.\tA\tATT\t30\tPASS\t.\tGT\t0/1\n')
    (work / 'base.bed').write_text('chrT\t10\t20\n')
    (work / 'old.bed').write_text('chrT\t19\t21\n')
    run('generate', [str(tools / 'gvcf2bed'), 'truth.vcf', '-r', 'ref.fa', '-o', 'fixed.bed'])
    if (work / 'fixed.bed').read_text() != 'chrT\t19\t22\tcall,insertion:PASS\n':
        raise ValueError('The union right endpoint for adjacent insertions must be 22')
    score_header = header + (
        '##INFO=<ID=BS,Number=1,Type=Integer,Description="Superlocus">\n'
        '##INFO=<ID=type,Number=1,Type=String,Description="Decision">\n'
        '##INFO=<ID=kind,Number=1,Type=String,Description="Match kind">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n')
    rows = [(11, 'C', 'TP', '0/1', '0/1'), (19, 'C', 'TP', '0/1', '0/1'),
            (20, 'ATT', 'TP', '0/1', '0/1'), (21, 'ATT', 'FP', '0/1', './.'),
            (22, 'C', 'FP', './.', '0/1')]
    (work / 'score.vcf').write_text(score_header + ''.join(
        f'chrT\t{pos}\t.\tA\t{alt}\t30\tPASS\tBS=11;type={typ};kind=almismatch\tGT\t{truth}\t{query}\n'
        for pos, alt, typ, truth, query in rows))
    pysam.tabix_compress(str(work / 'score.vcf'), str(work / 'score.vcf.gz'))
    pysam.tabix_index(str(work / 'score.vcf.gz'), preset='vcf')
    for version in ('old', 'fixed'):
        run(version, [str(tools / 'quantify'), 'score.vcf.gz', '--type', 'xcmp',
                      '-r', 'ref.fa', '-R', 'CONF:base.bed', '-R', f'CONF_VARS:{version}.bed',
                      '-o', f'{version}.tsv', '-v', f'{version}.vcf.gz', '--threads', '1'])
        pysam.tabix_index(str(work / f'{version}.vcf.gz'), preset='vcf')
        check_records(work / f'{version}.vcf.gz', fixed=version == 'fixed')
        with (work / f'{version}.tsv').open() as stream:
            table_rows = list(csv.DictReader(stream, delimiter='\t'))
        check_counts(table_rows, fixed=version == 'fixed')
    try:
        check_records(work / 'old.vcf.gz', fixed=True)
    except ValueError as exc:
        results.append({'name': 'legacy-output-rejected-by-fixed-contract', 'detected': True, 'reason': str(exc)})
    else:
        raise RuntimeError('Negative control did not detect the legacy boundary error')
    for field in ('BD', 'Regions'):
        negative_path = work / f'negative-{field}.vcf'
        with pysam.VariantFile(str(work / 'fixed.vcf.gz')) as source:
            with pysam.VariantFile(str(negative_path), 'w', header=source.header) as target:
                for record in source:
                    if record.pos == 22:
                        if field == 'BD':
                            record.samples['QUERY']['BD'] = 'UNK'
                        else:
                            record.info['Regions'] = ('TS_boundary',)
                    target.write(record)
        try:
            check_records(negative_path, fixed=True)
        except ValueError as exc:
            results.append({'name': field + '-mutation-rejected', 'detected': True, 'reason': str(exc)})
        else:
            raise RuntimeError(f'Negative control did not detect perturbation of {field}')
    mutated = [dict(row) for row in table_rows]
    for row in mutated:
        if row['Subset'] == row['Subtype'] == row['Genotype'] == row['QQ'] == '*' and row['Filter'] == 'ALL' and row['Type'] == 'SNP':
            row['TRUTH.TOTAL'] = '3'
    try:
        check_counts(mutated, fixed=True)
    except ValueError as exc:
        results.append({'name': 'count-mutation-rejected', 'detected': True, 'reason': str(exc)})
    else:
        raise RuntimeError('Negative control did not detect perturbed counts')
    (work / 'results.json').write_text(json.dumps({'contract': 'CONF-union-v1', 'results': results}, indent=2) + '\n')
    print('PASS: CONF-union-v1, five records, both insertion flanks, FN/FP, BS propagation, rejection of legacy output')


if __name__ == '__main__':
    main()
