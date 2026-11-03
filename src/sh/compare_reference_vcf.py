#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare VCFs strictly against the actual reference dictionary, preserving all other golden content."""
import argparse
from pathlib import Path


def reference_name(reference: Path) -> str:
    """Derive the expected reference identifier from the input filename."""
    name: str = reference.name
    for suffix in ('.fasta.gz', '.fa.gz', '.fasta', '.fa'):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def compare(actual: list[str], expected: list[str], reference: Path) -> None:
    fields: list[list[str]] = [line.split('\t') for line in
        Path(str(reference) + '.fai').read_text().splitlines()]
    contigs: list[str] = [f'##contig=<ID={row[0]},length={int(row[1])}>' for row in fields]
    if [line for line in actual if line.startswith('##contig=')] != contigs:
        raise AssertionError('Reference dictionary differs')
    if [line for line in actual if line.startswith('##reference=')] != [
            '##reference=' + reference_name(reference)]:
        raise AssertionError('Reference identifier differs')
    retained_actual: list[str] = [line for line in actual if not line.startswith(('##reference=', '##contig='))]
    retained_expected: list[str] = [line for line in expected if not line.startswith(('##reference=', '##contig='))]
    if retained_actual != retained_expected:
        raise AssertionError('Non-reference metadata or records differ')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('actual', type=Path)
    parser.add_argument('expected', type=Path)
    parser.add_argument('reference', type=Path)
    args = parser.parse_args()
    compare(args.actual.read_text().splitlines(), args.expected.read_text().splitlines(), args.reference)


if __name__ == '__main__':
    main()
