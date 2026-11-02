#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate original AD capacity, complete cross-chromosome reading, and failure on parse errors."""
import argparse
import json
import os
from pathlib import Path
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    work: Path = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    tools: Path = args.bin_dir.resolve()
    results: list[dict] = []

    def run(name: str, command: list[str], expected: int = 0) -> None:
        process = subprocess.run(command, cwd=work, capture_output=True, text=True, check=False)
        (work / f'{name}.stdout').write_text(process.stdout)
        (work / f'{name}.stderr').write_text(process.stderr)
        results.append({'name': name, 'command': command, 'exit_code': process.returncode})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        if expected == 0:
            if not (process.returncode == 0):
                raise AssertionError((name, process.stderr))
        else:
            if not (process.returncode != 0):
                raise AssertionError(name)
            if not ('Failed to read' in process.stderr):
                raise AssertionError(process.stderr)

    (work / 'ref.fa').write_text('>chr2\nACGTACGTACGT\n>chr10\nACGTACGTACGT\n')
    (work / 'conf.bed').write_text('chr2\t0\t12\nchr10\t0\t12\n')
    header: str = ('##fileformat=VCFv4.2\n##contig=<ID=chr10,length=12>\n'
                   '##contig=<ID=chr2,length=12>\n'
                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                   '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
    row: str = 'chr2\t2\t.\tC\tT\t30\tPASS\t.\tGT\t0/1\n'
    (work / 'input.vcf').write_text(header + row + row.replace('chr2', 'chr10'))
    (work / 'broken.vcf').write_text(header + row + row.replace('\t2\t', '\tbad\t', 1))
    run('faidx', ['samtools', 'faidx', 'ref.fa'])
    run('compress', ['bcftools', 'view', '-Oz', '-o', 'input.vcf.gz', 'input.vcf'])
    run('index', ['bcftools', 'index', '-t', 'input.vcf.gz'])
    expected: str = 'chr2\t1\t2\tcall:PASS\nchr10\t1\t2\tcall:PASS\n'
    for label, filename in [('plain', 'input.vcf'), ('indexed', 'input.vcf.gz')]:
        run(label, [str(tools / 'gvcf2bed'), filename, '-r', 'ref.fa', '-T', 'conf.bed', '-o', f'{label}.bed'])
        if not ((work / f'{label}.bed').read_text() == expected):
            raise AssertionError(label)
    (work / 'overlapping-insertions.vcf').write_text(
        header + 'chr2\t2\t.\tC\tCTT\t30\tPASS\t.\tGT\t0/1\n'
        + 'chr2\t3\t.\tG\tGTT\t30\tPASS\t.\tGT\t0/1\n')
    for label, options in [('overlap-plain', []), ('overlap-target', ['-T', 'conf.bed'])]:
        run(label, [str(tools / 'gvcf2bed'), 'overlapping-insertions.vcf', '-r', 'ref.fa',
                    *options, '-o', f'{label}.bed'])
        observed: str = (work / f'{label}.bed').read_text()
        expected_union: str = 'chr2\t1\t4\tcall,insertion:PASS\n'
        if observed != expected_union:
            raise ValueError(f'{label}: expected {expected_union!r}, observed {observed!r}')
    run('malformed', [str(tools / 'gvcf2bed'), 'broken.vcf', '-r', 'ref.fa', '-T', 'conf.bed', '-o', 'broken.bed'], expected=1)
    (work / 'ad.fa').write_text('>chr1\nA' + 'C' * 63 + '\n')
    ad_header: str = header.replace('##contig=<ID=chr10,length=12>\n##contig=<ID=chr2,length=12>', '##contig=<ID=chr1,length=64>')
    ad_header = ad_header.replace('#CHROM', '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n#CHROM')
    (work / 'ad.vcf').write_text(ad_header + 'chr1\t2\t.\tC\tCTTT,CTTT\t30\tPASS\t.\tGT:AD\t2/1:0,8,13\n')
    run('ad-faidx', ['samtools', 'faidx', 'ad.fa'])
    run('ad-compress', ['bcftools', 'view', '-Oz', '-o', 'ad.vcf.gz', 'ad.vcf'])
    run('ad-index', ['bcftools', 'index', 'ad.vcf.gz'])
    run('ad-preprocess', [str(tools / 'preprocess'), 'ad.vcf.gz:*', '-o', 'ad-output.vcf.gz', '-V', '1', '-L', '1', '-r', 'ad.fa'])
    run('ad-records', ['bcftools', 'view', '-H', 'ad-output.vcf.gz'])
    if not ((work / 'ad-records.stdout').read_text().strip()):
        raise AssertionError("Regression condition failed: (work / 'ad-records.stdout').read_text().strip()")
    launcher: Path = tools.parent / 'libexec/rtg-tools-install/rtg'
    fake: Path = work / 'fake-java/bin'
    fake.mkdir(parents=True)
    (fake / 'java').write_text('#!/bin/sh\nprintf "%s\\n" "$@"\n')
    (fake / 'java').chmod(0o755)
    env: dict[str, str] = dict(os.environ, JAVA_HOME=str(fake.parent), RTG_MEM='256M')
    process = subprocess.run([str(launcher), 'version'], env=env, capture_output=True, text=True, check=True)
    if not ('-Xmx256M' in process.stdout.splitlines()):
        raise AssertionError(process.stdout)
    if not ('-Dtalkback=false' in process.stdout.splitlines()):
        raise AssertionError(process.stdout)
    overlays: list[Path] = list((launcher.parent / 'lib').glob('htsjdk-*-rtg-overlay.jar'))
    if not (len(overlays) == 1):
        raise AssertionError(overlays)
    if not (str(overlays[0]) in process.stdout):
        raise AssertionError(process.stdout)
    (work / 'rtg-memory.stdout').write_text(process.stdout)
    results.append({'name': 'rtg-memory-argument', 'exit_code': process.returncode})
    (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
    print('PASS: AD capacity, complete plain/indexed VCF reading, malformed-input failure, RTG_MEM arguments')


if __name__ == '__main__':
    main()
