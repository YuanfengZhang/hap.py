#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Preserve legacy strict failures and independently run the versioned HTSJDK6 validation contract."""
import argparse
import json
from pathlib import Path
import re
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--java-home', type=Path, required=True)
    parser.add_argument('--classpath', required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    root: Path = Path(__file__).resolve().parent
    work: Path = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    classes: Path = work / 'classes'
    classes.mkdir()
    results: list[dict] = []

    def run(name: str, command: list[str]) -> subprocess.CompletedProcess:
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        (work / (name + '.log')).write_text(process.stdout + process.stderr)
        results.append({'name': name, 'command': command, 'exit_code': process.returncode})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        return process

    java: Path = args.java_home.resolve() / 'bin'
    compile_result = run('compile-contract', [str(java / 'javac'), '-cp', args.classpath,
        '-d', str(classes), str(root / 'patches/Htsjdk6ValidationContractTest.java'),
        str(root / 'patches/WriterDictionaryTest.java')])
    if not (compile_result.returncode == 0):
        raise AssertionError(compile_result.stderr)
    command: list[str] = [str(java / 'java'), '-Xmx8g', '-Drtg.developer=true',
        '-Dtalkback=false', '-Drtg.version.check-latest=false', '-Djava.io.tmpdir=' + str(work),
        '-cp', str(classes) + ':' + args.classpath, 'org.junit.runner.JUnitCore']
    strict = run('legacy-strict', command + (root / 'patches/rtg-regression.classes').read_text().splitlines())
    failures: list[str] = re.findall(r'^\d+\) (\S+)', strict.stdout + strict.stderr, re.MULTILINE)
    expected: list[str] = [
        'testBad2First5(com.rtg.sam.SamFileAndRecordTest)',
        'testBad2First5(com.rtg.sam.SkipInvalidRecordsIteratorTest)',
    ]
    if not (strict.returncode == 1 and sorted(failures) == sorted(expected)):
        raise AssertionError(failures)
    strict_text: str = strict.stdout + strict.stderr
    if not (re.search(r'^Tests run: 1001,  Failures: 2$', strict_text, re.MULTILINE)):
        raise AssertionError(strict_text)
    blocks: list[str] = re.split(r'^\d+\) ', strict_text, flags=re.MULTILINE)[1:]
    if not (len(blocks) == 2):
        raise AssertionError(blocks)
    if not (sorted(block.splitlines()[0] for block in blocks) == sorted(expected)):
        raise AssertionError(blocks)
    for block in blocks:
        if not ("java.lang.AssertionError: 'ERROR::INVALID_FLAG_NOT_PRIM_ALIGNMENT:" in block):
            raise AssertionError(block)
        if not ('was not contained in:SAM record is invalid in file ' in block):
            raise AssertionError(block)
        if not ('WARNING::INVALID_FLAG_NOT_PRIM_ALIGNMENT:' in block):
            raise AssertionError(block)
        if not ('WARNING::INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT:' in block):
            raise AssertionError(block)
        if not ('At data line 1' in block):
            raise AssertionError(block)
        if not ('INVALID_UNMAPPED_MAPPING_QUALITY' not in block):
            raise AssertionError(block)
    contract = run('versioned-contract', command + ['com.rtg.sam.Htsjdk6ValidationContractTest',
                                                  'com.rtg.sam.WriterDictionaryTest'])
    if not (contract.returncode == 0):
        raise AssertionError(contract.stdout + contract.stderr)
    (work / 'acceptance.json').write_text(json.dumps({
        'legacy_strict_exit_code': strict.returncode,
        'legacy_strict_pass': False,
        'expected_upstream_contract_differences': failures,
        'versioned_contract_exit_code': contract.returncode,
        'scope': 'Validate only the listed upstream semantic differences; do not claim the legacy strict suite passes.',
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
