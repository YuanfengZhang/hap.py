#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate the distance-v1 neighborhood contract without treating distance TPs as haplotype equivalence."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import pysam


def decisions(path: Path) -> list[tuple[str, str]]:
    with pysam.VariantFile(str(path)) as stream:
        return [tuple(record.samples[name].get('BD') for name in ('TRUTH', 'QUERY')) for record in stream]


def check_identity(source: Path, output: Path) -> None:
    def rows(path: Path) -> tuple[tuple[str, ...], list[tuple]]:
        with pysam.VariantFile(str(path)) as stream:
            names = tuple(stream.header.samples)
            records = [(record.contig, record.pos, record.ref, record.alts, tuple(record.filter),
                        tuple((name, record.samples[name].get('GT'), record.samples[name].phased) for name in names))
                       for record in stream]
        return names, records
    if rows(source) != rows(output):
        raise ValueError('distance output changed record identity, order, FILTER, or GT')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    tools, work = args.bin_dir.resolve(), args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    reference = work / 'ref.fa'
    reference.write_text('>chrT\n' + 'A' * 4096 + '\n>chrU\n' + 'A' * 4096 + '\n')
    pysam.faidx(str(reference))
    header = ('##fileformat=VCFv4.2\n##contig=<ID=chrT,length=4096>\n##contig=<ID=chrU,length=4096>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n')
    cases = {
        'inclusive-distance': ([('chrT', 100, 'C', '0/1', './.'), ('chrT', 130, 'G', './.', '0/1')], [('TP', '.'), ('.', 'TP')]),
        'outside-distance': ([('chrT', 100, 'C', '0/1', './.'), ('chrT', 131, 'G', './.', '0/1')], [('FN', '.'), ('.', 'FP')]),
        'alt-span': ([('chrT', 100, 'A' + 'T' * 40, '0/1', './.'), ('chrT', 170, 'C', './.', '0/1')], [('TP', '.'), ('.', 'TP')]),
        'alt-span-outside': ([('chrT', 100, 'A' + 'T' * 40, '0/1', './.'), ('chrT', 171, 'C', './.', '0/1')], [('FN', '.'), ('.', 'FP')]),
        'chromosome-isolation': ([('chrT', 100, 'C', '0/1', './.'), ('chrU', 100, 'G', './.', '0/1')], [('FN', '.'), ('.', 'FP')]),
        'many-to-one': ([('chrT', 100, 'C', '0/1', './.'), ('chrT', 101, 'G', '0/1', './.'), ('chrT', 130, 'T', './.', '0/1')], [('TP', '.'), ('TP', '.'), ('.', 'TP')]),
    }
    results: list[dict] = []
    for name, (rows, expected) in cases.items():
        source, output = work / f'{name}.vcf', work / f'{name}.out.vcf'
        source.write_text(header + ''.join(f'{chrom}\t{pos}\t.\tA\t{alt}\t30\tPASS\t.\tGT\t{truth}\t{query}\n' for chrom, pos, alt, truth, query in rows))
        with source.open() as stream:
            process = subprocess.run([str(tools / 'scmp'), '-M', 'distance', '-', '-r', str(reference),
                                      '--distance-maxdist', '30', '--threads', '1', '-o', str(output)],
                                     stdin=stream, capture_output=True, text=True, check=False)
        (work / f'{name}.log').write_text(process.stdout + process.stderr)
        if process.returncode == 0:
            check_identity(source, output)
        observed = decisions(output) if process.returncode == 0 else []
        passed = process.returncode == 0 and observed == expected
        results.append({'case': name, 'compute_exit_code': process.returncode, 'expected': expected,
                        'observed': observed, 'passed': passed})
        (work / 'results.json').write_text(json.dumps({'contract': 'distance-v1', 'results': results}, indent=2) + '\n')
    if not all(row['passed'] for row in results):
        raise SystemExit(1)
    if decisions(work / 'outside-distance.out.vcf') == cases['inclusive-distance'][1]:
        raise RuntimeError('Out-of-threshold negative control did not fail')
    print('PASS: distance-v1, threshold boundaries, ALT spans, chromosome isolation, many-to-one matching; not proof of haplotype equivalence')


if __name__ == '__main__':
    main()
