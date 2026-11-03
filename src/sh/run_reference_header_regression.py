#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate reference-dictionary migration while preserving all other headers and variant records verbatim."""
import argparse
import json
from pathlib import Path
import subprocess


def reference_name(reference: Path) -> str:
    """Derive the expected reference identifier from the input filename."""
    name: str = reference.name
    for suffix in ('.fasta.gz', '.fa.gz', '.fasta', '.fa'):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin-dir', type=Path, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    root: Path = Path(__file__).resolve().parents[2]
    work: Path = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    cases: list[tuple[str, str, str, list[str]]] = [
        ('merge', 'src/data/microhg19.fa', 'src/data/expected_merge.vcf',
         [str(root / 'src/data/merge1.vcf.gz') + ':NA12877',
          str(root / 'src/data/merge2.vcf.gz') + ':NA12878',
          '--trimalleles=1', '--merge-by-location=1', '--unique-alleles=1']),
        ('callsonly', 'example/chr21.fa', 'example/callsonly/expected_callsonly.vcf',
         [str(root / 'example/callsonly/call_merge.vcf.gz') + ':*',
          '--process-full=1', '--process-formats=1']),
    ]
    results: list[dict] = []
    for name, reference, golden, options in cases:
        output: Path = work / (name + '.vcf')
        command: list[str] = [str(args.bin_dir.resolve() / 'multimerge'), *options,
                              '-r', str(root / reference), '-o', str(output)]
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        (work / (name + '.stderr')).write_text(process.stderr)
        result: dict = {'name': name, 'command': command, 'exit_code': process.returncode}
        results.append(result)
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        if not (process.returncode == 0):
            raise AssertionError(process.stderr)
        actual: list[str] = output.read_text().splitlines()
        expected: list[str] = (root / golden).read_text().splitlines()
        fai: list[tuple[str, int]] = [(fields[0], int(fields[1])) for fields in
            (line.split('\t') for line in Path(str(root / reference) + '.fai').read_text().splitlines())]
        result['reference_dictionary_exact'] = [line for line in actual if line.startswith('##contig=')] == [
            f'##contig=<ID={name},length={length}>' for name, length in fai]
        result['reference_name_exact'] = [line for line in actual if line.startswith('##reference=')] == [
            '##reference=' + reference_name(root / reference)]
        retained_actual = [line for line in actual if not line.startswith(('##reference=', '##contig='))]
        retained_expected = [line for line in expected if not line.startswith(('##reference=', '##contig='))]
        result['other_metadata_and_records_exact'] = retained_actual == retained_expected
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        if not (all(result[key] for key in ('reference_dictionary_exact', 'reference_name_exact',
                                          'other_metadata_and_records_exact'))):
            raise AssertionError(result)
    print('PASS: actual reference dictionary, reference identifier, and all other metadata and variant records')


if __name__ == '__main__':
    main()
