#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify quarantine CLI output against independently filtered native inputs.

Run after test_native_engines.py, which stages the current Python sources.
Usage: python -I -B test_invalid_native.py NATIVE_BIN
"""
import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import pysam

ROOT = Path(__file__).resolve().parents[3]


def main() -> None:
    native = Path(sys.argv[1]).resolve()
    stage = ROOT / '.work/installed-integration'
    environment = dict(os.environ, PATH=str(native) + os.pathsep + os.environ['PATH'])
    header = ('##fileformat=VCFv4.2\n##contig=<ID=c,length=400>\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
              '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n')
    with tempfile.TemporaryDirectory(prefix='quarantine é ') as directory:
        root = Path(directory)
        reference = root / 'ref.fa'
        reference.write_text('>c\n' + 'ACGT' * 100 + '\n')
        pysam.faidx(str(reference))
        good = 'c\t2\t.\tC\tT\t30\tPASS\t.\tGT:AD\t0/1:8,2\n'
        bad = ('c\t10\t.\tC\tT\t30\tPASS\t.\tGT:AD\t0/2:8,2\n'
               'c\t14\t.\tC\tT\t30\tPASS\t.\tGT:AD\t0/1/1:8,2\n'
               'c\t18\t.\tC\tT\t30\tPASS\t.\tGT:AD\t0/1:8,2,3\n')
        bed = root / 'conf.bed'
        bed.write_text('c\t0\t8\nc\t20\t100\n')
        paths: dict[str, Path] = {}
        for label in ('truth', 'query'):
            for mode in ('clean', 'mixed'):
                path = root / f'{label}.{mode}.vcf.gz'
                with pysam.BGZFile(str(path), 'w') as out:
                    extra = good.replace('\t2\t', '\t22\t' if label == 'truth' else '\t26\t')
                    out.write((header + good + (bad if mode == 'mixed' else '') + extra).encode())
                pysam.tabix_index(str(path), preset='vcf')
                paths[label + mode] = path

        def run(command: list[str], success: bool = True) -> subprocess.CompletedProcess:
            result = subprocess.run(command, env=environment, capture_output=True, text=True)
            if (result.returncode == 0) != success:
                raise AssertionError((command, result.returncode, result.stdout, result.stderr))
            return result

        def records(path: Path) -> list[str]:
            with pysam.VariantFile(path) as source:
                return [str(record) for record in source]

        pre = [sys.executable, '-I', '-B', str(stage / 'bin/pre.py')]
        for mode in ('clean', 'mixed'):
            output = root / f'{mode}.pre.vcf.gz'
            run(pre + [str(paths['query' + mode]), str(output), '-r', str(reference),
                       '--threads', '2'] + (['--output-invalid'] if mode == 'mixed' else []))
        assert records(root / 'clean.pre.vcf.gz') == records(root / 'mixed.pre.vcf.gz')
        assert records(root / 'mixed.pre.vcf.gz.invalid.vcf.gz') == bad.splitlines(keepends=True)
        run(pre + [str(paths['querymixed']), str(root / 'strict.vcf.gz'), '-r', str(reference)], False)

        all_bad = root / 'all-bad.vcf'
        all_bad.write_text(header + bad)
        all_output = root / 'all-bad.pre.vcf.gz'
        failed = run(pre + [str(all_bad), str(all_output), '-r', str(reference),
                            '--output-invalid'], False)
        assert 'No valid records remain' in failed.stderr
        assert records(Path(str(all_output) + '.invalid.vcf.gz')) == bad.splitlines(keepends=True)
        failed = run(pre + [str(paths['querymixed']), str(root / 'mixed.pre.vcf.gz'),
                            '-r', str(reference), '--output-invalid'], False)
        assert 'Refusing to overwrite' in failed.stderr
        assert not list(Path(os.environ['TMPDIR']).glob('valid.*'))

        chr_ref = root / 'chr-ref.fa'
        chr_ref.write_text('>chr1\n' + 'ACGT' * 100 + '\n')
        pysam.faidx(str(chr_ref))
        bare = root / 'bare.vcf'
        bare.write_text((header + good + bad).replace('ID=c,', 'ID=1,').replace('\nc\t', '\n1\t'))
        chr_output = root / 'chr.pre.vcf.gz'
        run(pre + [str(bare), str(chr_output), '-r', str(chr_ref), '--output-invalid'])
        assert records(chr_output)[0].startswith('chr1\t')

        for engine in ('xcmp', 'scmp-somatic', 'vcfeval'):
            for threads in (1, 2):
                prefixes: list[Path] = []
                for mode in ('clean', 'mixed'):
                    prefix = root / f'{engine}.{threads}.{mode}'
                    prefixes.append(prefix)
                    command = [sys.executable, '-I', '-B', str(stage / 'bin/hap.py'),
                               str(paths['truth' + mode]), str(paths['query' + mode]),
                               '-r', str(reference), '-f', str(bed), '-o', str(prefix), '--engine', engine,
                               '--threads', str(threads)]
                    if engine == 'vcfeval':
                        command.extend(['--engine-vcfeval-path',
                                        str(native.parent / 'libexec/rtg-tools-install/rtg')])
                    if mode == 'mixed':
                        command.append('--output-invalid')
                    run(command)
                    if mode == 'mixed':
                        for label in ('truth', 'query'):
                            assert records(Path(str(prefix) + f'.{label}.invalid.vcf.gz')) == bad.splitlines(keepends=True)
                            reasons = Path(str(prefix) + f'.{label}.invalid.vcf.gz.reasons.tsv').read_text()
                            assert all(reason in reasons for reason in ('GT_ALLELE', 'GT_PLOIDY', 'AD_CAPACITY'))
                clean, mixed = prefixes
                for table in root.glob(clean.name + '.*'):
                    suffix = table.name[len(clean.name):]
                    if not suffix.endswith(('.csv', '.csv.gz')):
                        continue
                    other = Path(str(mixed) + suffix)
                    if suffix.endswith('.gz'):
                        with gzip.open(table, 'rt') as left, gzip.open(other, 'rt') as right:
                            assert left.read() == right.read(), suffix
                    else:
                        assert table.read_text() == other.read_text(), suffix
                assert records(Path(str(clean) + '.vcf.gz')) == records(Path(str(mixed) + '.vcf.gz'))
                print(f'PASS {engine} threads={threads}: tables, ROC, VCF, quarantine', flush=True)
        print('PASS pre.py strict/quarantine and native-engine equivalence', flush=True)


if __name__ == '__main__':
    main()
