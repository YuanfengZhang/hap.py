#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check scmp allele matching against explicit reference-span expectations."""
import argparse
import json
import subprocess
from pathlib import Path

import pysam


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--bin-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=False)
    reference = args.work_dir / 'reference.fa'
    reference.write_text('>chrT\n' + 'AC' * 50 + 'GGAAAGAAAAAGAAAGAAAGAAA' + 'CT' * 100 + '\n')
    pysam.faidx(str(reference))
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chrT,length=321>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n')
    cases = {
        'padded_deletions': ([
            ('GGAAAGAAAAA', 'G', './1', './1'),
            ('GGAAAGAAAAAGAAAGAAAGAAA', 'GGAAA,G', './2', './1'),
            ('GGAAAGAAAAAGAAAGAAA', 'G', './1', './.')], {'TP': 4, 'FN': 1, 'FP': 0}),
        'different_reference_span': ([
            ('GGAAAGAAAAA', 'C', './1', './.'),
            ('G', 'C', './.', './1')], {'TP': 0, 'FN': 1, 'FP': 1}),
        'equal_substitution': ([('G', 'C', './1', './1')], {'TP': 2, 'FN': 0, 'FP': 0}),
        'padded_insertion': ([
            ('G', 'GTT', './1', './.'),
            ('GG', 'GTTG', './.', './1')], {'TP': 2, 'FN': 0, 'FP': 0}),
        'different_insertions': ([
            ('G', 'GTT', './1', './.'),
            ('G', 'GCC', './.', './1')], {'TP': 0, 'FN': 1, 'FP': 1}),
    }
    expected_rows = {
        'padded_deletions': [('TP', 'TP'), ('FN', 'TP'), ('TP', '.')],
        'different_reference_span': [('FN', '.'), ('.', 'FP')],
        'equal_substitution': [('TP', 'TP')],
        'padded_insertion': [('TP', '.'), ('.', 'TP')],
        'different_insertions': [('FN', '.'), ('.', 'FP')],
    }
    results: list[dict] = []
    for name, (records, expected) in cases.items():
        source = args.work_dir / f'{name}.vcf'
        source.write_text(header + ''.join(
            f'chrT\t101\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{truth}\t{query}\n'
            for ref, alt, truth, query in records))
        output = args.work_dir / f'{name}.scored.vcf'
        with source.open() as stream, (args.work_dir / f'{name}.log').open('w') as log:
            process = subprocess.run([str(args.bin_dir / 'scmp'), '-M', 'alleles', '-',
                                      '-r', str(reference), '--threads', '1', '-o', str(output)],
                                     stdin=stream, stdout=log, stderr=subprocess.STDOUT, check=False)
        observed = dict.fromkeys(expected, 0)
        observed_rows: list[tuple] = []
        if process.returncode == 0:
            with pysam.VariantFile(str(output)) as vcf:
                for record in vcf:
                    observed_rows.append(tuple(record.samples[n].get('BD') for n in ('TRUTH', 'QUERY')))
                    for sample in record.samples.values():
                        decision = sample.get('BD')
                        if decision in observed:
                            observed[decision] += 1
        results.append({'case': name, 'compute_exit_code': process.returncode,
                        'expected': expected, 'observed': observed,
                        'expected_rows': expected_rows[name], 'observed_rows': observed_rows,
                        'passed': process.returncode == 0 and observed == expected
                        and observed_rows == expected_rows[name]})
    for name, options in (('default_mode', []), ('singular_alias', ['-M', 'allele'])):
        source = args.work_dir / 'equal_substitution.vcf'
        output = args.work_dir / f'{name}.scored.vcf'
        with source.open() as stream, (args.work_dir / f'{name}.log').open('w') as log:
            process = subprocess.run([str(args.bin_dir / 'scmp'), *options, '-',
                                      '-r', str(reference), '--threads', '1', '-o', str(output)],
                                     stdin=stream, stdout=log, stderr=subprocess.STDOUT, check=False)
        observed_rows: list[tuple] = []
        if process.returncode == 0:
            with pysam.VariantFile(str(output)) as vcf:
                observed_rows = [tuple(record.samples[n].get('BD') for n in ('TRUTH', 'QUERY'))
                                 for record in vcf]
        results.append({'case': name, 'compute_exit_code': process.returncode,
                        'expected_rows': expected_rows['equal_substitution'],
                        'observed_rows': observed_rows,
                        'passed': process.returncode == 0 and
                        observed_rows == expected_rows['equal_substitution']})
    for mode in ('alleles', 'enumerate'):
        for index, alt in enumerate(('<DEL>', '<INS>', '*', 'G]chrT:120]')):
            source = args.work_dir / f'unsupported-{mode}-{index}.vcf'
            source.write_text(header + f'chrT\t101\t.\tG\t{alt}\t.\tPASS\t.\tGT\t./1\t./1\n')
            with source.open() as stream:
                process = subprocess.run([str(args.bin_dir / 'scmp'), '-M', mode, '-',
                                          '-r', str(reference), '--threads', '1', '-o',
                                          str(args.work_dir / f'unsupported-{mode}-{index}.out.vcf')],
                                         stdin=stream, capture_output=True, text=True, check=False)
            (args.work_dir / f'unsupported-{mode}-{index}.log').write_text(process.stdout + process.stderr)
            results.append({'case': f'unsupported-{mode}-{alt}', 'compute_exit_code': process.returncode,
                            'passed': process.returncode != 0 and
                            'unsupported non-sequence ALT' in process.stderr})
    (args.work_dir / 'result.json').write_text(json.dumps(results, indent=2) + '\n')
    failures = [result['case'] for result in results if not result['passed']]
    if failures:
        raise ValueError('Reference-span regressions failed: ' + ', '.join(failures))


if __name__ == '__main__':
    main()
