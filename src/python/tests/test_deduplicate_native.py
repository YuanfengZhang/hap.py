#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exercise real normalization collisions and duplicate compatibility switch."""
from pathlib import Path
import os
import gzip
import subprocess
import sys
import tempfile

import pysam

ROOT = Path(__file__).resolve().parents[3]


def main() -> None:
    native = Path(sys.argv[1]).resolve()
    stage = ROOT / '.work/installed-integration'
    env = dict(os.environ, PATH=str(native) + os.pathsep + os.environ['PATH'])
    with tempfile.TemporaryDirectory(prefix='dedup native é ') as folder:
        root = Path(folder)
        reference = root / 'ref.fa'
        reference.write_text('>c\nC' + 'A' * 100 + 'G\n')
        pysam.faidx(str(reference))
        header = ('##fileformat=VCFv4.2\n##contig=<ID=c,length=102>\n'
                  '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n')
        rows = 'c\t5\t.\tA\tAA\t30\tPASS\t.\tGT\t0/1\nc\t8\t.\tA\tAA\t30\tPASS\t.\tGT\t0/1\n'
        source = root / 'input.vcf.gz'
        with pysam.BGZFile(str(source), 'w') as output:
            output.write((header + rows).encode())
        pysam.tabix_index(str(source), preset='vcf')

        def run(arguments: list[str]) -> None:
            p = subprocess.run(arguments, env=env, capture_output=True, text=True)
            if p.returncode:
                raise AssertionError((arguments, p.stdout, p.stderr))

        # Independent normalization, without bcftools allele-only deduplication.
        normalized = root / 'normalized.vcf.gz'
        run(['bcftools', 'norm', '-f', str(reference), '-Oz', '-o', str(normalized), str(source)])
        with pysam.VariantFile(normalized) as records:
            text = [str(record) for record in records]
        assert len(text) == 2 and text[0] == text[1], text
        print('PASS: distinct repetitive-sequence input positions normalize to identical records', flush=True)
        # Feed normalized collision through final CLI output, not merely a mock.
        pysam.tabix_index(str(normalized), preset='vcf')
        for threads in (1, 2):
            for bcf in (False, True):
                for keep in (False, True):
                    output = root / f'out-{threads}-{bcf}-{keep}.vcf.gz'
                    command = [sys.executable, '-I', '-B', str(stage / 'bin/pre.py'),
                               str(normalized), str(output), '-r', str(reference), '--threads', str(threads),
                               '--no-leftshift', '--no-decompose', '--gender', 'none']
                    if bcf:
                        command.append('--bcf')
                    if keep:
                        command.append('--keep-duplicates')
                    run(command)
                    actual = Path(str(output) + '.bcf') if bcf else output
                    with pysam.VariantFile(actual) as records:
                        assert len(list(records)) == (2 if keep else 1)
        clean = root / 'clean.vcf.gz'
        with pysam.VariantFile(normalized) as source_records:
            with pysam.VariantFile(clean, 'wz', header=source_records.header) as dest:
                dest.write(next(source_records))
        pysam.tabix_index(str(clean), preset='vcf')
        for engine in ('xcmp', 'scmp-somatic', 'vcfeval'):
            for threads in (1, 2):
                prefixes: list[Path] = []
                for label, source_path in (('clean', clean), ('duplicate', normalized)):
                    prefix = root / f'{engine}-{threads}-{label}'
                    prefixes.append(prefix)
                    command = [sys.executable, '-I', '-B', str(stage / 'bin/hap.py'),
                               str(source_path), str(source_path), '-r', str(reference), '-o', str(prefix),
                               '--engine', engine, '--threads', str(threads), '--gender', 'none',
                               '--no-leftshift', '--no-decompose']
                    if engine == 'vcfeval':
                        command.extend(['--engine-vcfeval-path', str(native.parent / 'libexec/rtg-tools-install/rtg')])
                    run(command)
                left, right = prefixes
                for table in root.glob(left.name + '.*'):
                    suffix = table.name[len(left.name):]
                    if suffix.endswith('.csv'):
                        assert table.read_text() == Path(str(right) + suffix).read_text(), suffix
                    elif suffix.endswith('.csv.gz'):
                        with gzip.open(table, 'rt') as a, gzip.open(str(right) + suffix, 'rt') as b:
                            assert a.read() == b.read(), suffix
                with pysam.VariantFile(str(left) + '.vcf.gz') as a, pysam.VariantFile(str(right) + '.vcf.gz') as b:
                    assert [str(r) for r in a] == [str(r) for r in b]
                print(f'PASS {engine} threads={threads}: duplicate vs independent clean baseline', flush=True)
        # The norm path must retain same-allele records with different GT.
        with pysam.BGZFile(str(source), 'w') as output:
            output.write((header + rows.replace('c\t8', 'c\t5').replace('GT\t0/1\n', 'GT\t1/1\n', 1)).encode())
        pysam.tabix_index(str(source), preset='vcf', force=True)
        output = root / 'norm-different-gt.vcf.gz'
        run([sys.executable, '-I', '-B', str(stage / 'bin/pre.py'), str(source), str(output),
             '-r', str(reference), '--no-leftshift', '--no-decompose',
             '--bcftools-norm', '--gender', 'none'])
        with pysam.VariantFile(output) as records:
            assert len(list(records)) == 2
        print('PASS: threads/formats/keep switch and normalization retains distinct GT', flush=True)


if __name__ == '__main__':
    main()
