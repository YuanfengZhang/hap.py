#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exercise distance matching invariance across both historical block limits.

Usage: python -I test_scmp_distance_blocks.py SCMP_EXECUTABLE
"""
from pathlib import Path
import json
import os
import subprocess
import sys
import tempfile

import pysam


def main() -> None:
    executable = Path(sys.argv[1]).resolve()
    root = Path(__file__).resolve().parents[3]
    work = root / '.work'
    work.mkdir(exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='scmp-distance-regression-', dir=work))
    reference = directory / 'reference.fa'
    reference.write_text('>test\n' + 'C' * 50000 + '\n>other\n' + 'C' * 1000 + '\n')
    pysam.faidx(str(reference))
    header = ('##fileformat=VCFv4.2\n##contig=<ID=test,length=50000>\n'
              '##contig=<ID=other,length=1000>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n')
    results: list[dict] = []
    cases = [(f'alleles-{n}', n, '0/1') for n in (1022, 1023, 1024, 1025)]
    cases += [(f'records-{n}', n, '0/0') for n in (19999, 20000, 20001)]
    for name, count, genotype in cases:
        source = directory / f'{name}.vcf'
        with source.open('w') as stream:
            stream.write(header)
            for index in range(count):
                stream.write(f'test\t{index + 1}\t.\tC\tA\t40\tPASS\t.\tGT\t{genotype}\t./.\n')
            stream.write('test\t30000\t.\tCCCCCCCCCCCCC\tC\t40\tPASS\t.\tGT\t./.\t0/1\n')
            stream.write('test\t30004\t.\tC\tA\t40\tPASS\t.\tGT\t0/1\t./.\n')
            stream.write('other\t10\t.\tC\tA\t40\tPASS\t.\tGT\t0/1\t0/1\n')
        for threads in (1, 4):
            output = directory / f'{name}-{threads}.vcf'
            command = [str(executable), '-M', 'distance', '-', '-r', str(reference),
                       '--threads', str(threads), '--distance-maxdist', '30', '-o', str(output)]
            with source.open('rb') as input_stream, (directory / f'{name}-{threads}.log').open('w') as log:
                run = subprocess.run(command, stdin=input_stream, stdout=log, stderr=subprocess.STDOUT,
                                     env=os.environ.copy(), check=False)
            if not (run.returncode == 0):
                raise AssertionError((command, run.returncode))
            with pysam.VariantFile(str(output)) as variants:
                decisions = {(record.contig, record.pos): tuple(sample.get('BD') for sample in record.samples.values())
                             for record in variants if record.pos in (30000, 30004, 10)}
            if not (decisions[('test', 30000)][1] == 'TP'):
                raise AssertionError((name, threads, decisions))
            if not (decisions[('test', 30004)][0] == 'TP'):
                raise AssertionError((name, threads, decisions))
            if not (decisions[('other', 10)] == ('TP', 'TP')):
                raise AssertionError((name, threads, decisions))
            results.append({'case': name, 'threads': threads, 'exit_code': run.returncode})
    edge_cases = {
        'filtered-chromosome': (
            'test\t400\t.\tC\tA\t40\tPASS\t.\tGT\t0/1\t0/1\n'
            'other\t1\t.\tC\tA\t40\tLow\t.\tGT\t0/1\t0/1\n'
            'other\t2\t.\tC\tA\t40\tPASS\t.\tGT\t0/1\t0/1\n',
            ['--apply-filters', '1']),
        'symbolic-del': (
            'test\t101\t.\tC\tA\t40\tPASS\t.\tGT\t0/1\t./.\n'
            'test\t152\t.\tC\t<DEL>\t40\tPASS\t.\tGT\t./.\t0/1\n', []),
        'symbolic-del-same-position': (
            'test\t101\t.\tC\tA\t40\tPASS\t.\tGT\t0/1\t./.\n'
            'test\t152\t.\tC\tG\t40\tPASS\t.\tGT\t0/0\t0/0\n'
            'test\t152\t.\tC\t<DEL>\t40\tPASS\t.\tGT\t./.\t0/1\n', []),
    }
    for name, (rows, extra) in edge_cases.items():
        source = directory / f'{name}.vcf'
        source.write_text(header.replace('##FORMAT=', '##FILTER=<ID=Low,Description="Low">\n##FORMAT=', 1) + rows)
        for threads in (1, 4):
            output = directory / f'{name}-{threads}.vcf'
            command = [str(executable), '-M', 'distance', '-', '-r', str(reference),
                       '--threads', str(threads), '--distance-maxdist', '50', '-o', str(output), *extra]
            with source.open('rb') as input_stream, (directory / f'{name}-{threads}.log').open('w') as log:
                run = subprocess.run(command, stdin=input_stream, stdout=log, stderr=subprocess.STDOUT, check=False)
            if not (run.returncode == 0):
                raise AssertionError((command, run.returncode))
            with pysam.VariantFile(str(output)) as variants:
                records = list(variants)
            if name == 'filtered-chromosome':
                if not ([(r.contig, r.pos) for r in records] == [('test', 400), ('other', 2)]):
                    raise AssertionError("Regression condition failed: [(r.contig, r.pos) for r in records] == [('test', 400), ('other', 2)]")
                if not (all(tuple(s.get('BD') for s in r.samples.values()) == ('TP', 'TP') for r in records)):
                    raise AssertionError("Regression condition failed: all(tuple(s.get('BD') for s in r.samples.values()) == ('TP', 'TP') for r in records)")
            else:
                if not (records[0].samples['TRUTH']['BD'] == 'TP'):
                    raise AssertionError("Regression condition failed: records[0].samples['TRUTH']['BD'] == 'TP'")
                if not (records[-1].samples['QUERY']['BD'] == 'TP'):
                    raise AssertionError("Regression condition failed: records[-1].samples['QUERY']['BD'] == 'TP'")
            results.append({'case': name, 'threads': threads, 'exit_code': run.returncode})
    (directory / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    print(f'{len(results)} distance block invariance checks PASS: {directory}', flush=True)


if __name__ == '__main__':
    main()
