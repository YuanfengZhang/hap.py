#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The reference-dictionary contract must reject all unvalidated contig lines and data changes."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class ReferenceComparatorTest(unittest.TestCase):
    def test_optimized_reference_contract(self) -> None:
        result = subprocess.run(
            [sys.executable, '-I', '-O', str(Path(__file__).resolve()),
             'ReferenceComparatorTest.test_reference_contract'],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_reference_contract(self) -> None:
        root: Path = Path(__file__).resolve().parents[3]
        (root / '.work').mkdir(exist_ok=True)
        spec = importlib.util.spec_from_file_location('reference_comparator', root / 'src/sh/compare_reference_vcf.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory(dir=root / '.work', prefix='reference-contract-') as directory:
            reference: Path = Path(directory) / 'ref.fa'
            Path(str(reference) + '.fai').write_text('chrA\t10\t0\t10\t11\nchrB\t20\t0\t20\t21\n')
            actual: list[str] = ['##reference=ref', '##contig=<ID=chrA,length=10>',
                                 '##contig=<ID=chrB,length=20>', '#CHROM\tPOS', 'chrA\t1']
            expected: list[str] = ['##reference=hg19', '##contig=<ID=old,length=1>', '#CHROM\tPOS', 'chrA\t1']
            module.compare(actual, expected, reference)
            mutations: list[list[str]] = [
                ['##reference=wrong', *actual[1:]],
                actual + ['##contig=<ID=extra>'],
                actual + ['##contig=<ID=extra,length=1,assembly=x>'],
                actual + [actual[1]],
                [actual[0], actual[2], actual[1], *actual[3:]],
                [line.replace('length=10', 'length=11') for line in actual],
                [*actual[:-1], 'chrA\t2'],
                [*actual, '##INFO=<ID=extra>'],
            ]
            for index, changed in enumerate(mutations):
                with self.subTest(index=index), self.assertRaises(AssertionError):
                    module.compare(changed, expected, reference)


if __name__ == '__main__':
    unittest.main()
