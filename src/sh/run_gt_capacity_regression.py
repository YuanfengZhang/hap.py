#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify unsupported ploidy is rejected and valid multiallelic diploid genotypes are preserved."""
import argparse
import json
from pathlib import Path
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    root: Path = Path(__file__).resolve().parents[2]
    work: Path = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    reference: Path = root / 'src/data/chrQ.fa'
    original: str = (root / 'src/data/import_errors.vcf').read_text()
    results: list[dict] = []
    for name, gt, expected_success in [('triploid', '0/1/2', False), ('diploid', '1/2', True)]:
        source: Path = work / (name + '.vcf')
        source.write_text(original.replace('0/1/2', gt) + '\n')
        compressed: Path = work / (name + '.vcf.gz')
        subprocess.run(['bcftools', 'view', '-Oz', '-o', str(compressed), str(source)], check=True)
        subprocess.run(['bcftools', 'index', '-t', str(compressed)], check=True)
        output: Path = work / (name + '.output.vcf')
        command: list[str] = [str(args.bin_dir.resolve() / 'multimerge'), str(compressed) + ':NA12877',
                              '-r', str(reference), '-o', str(output)]
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        (work / (name + '.stderr')).write_text(process.stderr)
        rows: list[str] = [line for line in output.read_text().splitlines() if line and not line.startswith('#')] if output.exists() else []
        results.append({'name': name, 'command': command, 'exit_code': process.returncode,
                        'records': rows})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        if expected_success:
            if not (process.returncode == 0):
                raise AssertionError(process.stderr)
            if not (len(rows) == 1):
                raise AssertionError(rows)
            fields: list[str] = rows[0].split('\t')
            if not (fields[:5] == ['chrQ', '1', '.', 'AAAC', 'A,AA,AAA']):
                raise AssertionError(fields)
            if not ('IMPORT_FAIL' not in fields[7]):
                raise AssertionError(fields)
            if not (fields[9].split(':')[fields[8].split(':').index('GT')] == '1/2'):
                raise AssertionError(fields)
        else:
            if not (process.returncode != 0):
                raise AssertionError('Unsupported GT was accepted')
            if not ('GT ploidy exceeds supported capacity at chrQ:1' in process.stderr):
                raise AssertionError(process.stderr)
            if not (not rows):
                raise AssertionError(rows)
    print('PASS: triploid genotypes explicitly rejected; valid multiallelic diploid genotypes fully preserved')


if __name__ == '__main__':
    main()
