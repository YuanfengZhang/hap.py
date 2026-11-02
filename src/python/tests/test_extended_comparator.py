#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extended-table comparisons preserve exact counts, complete columns, and stratification contracts."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class ExtendedComparatorTest(unittest.TestCase):
    def test_numeric_and_schema_contract(self) -> None:
        root: Path = Path(__file__).resolve().parents[3]
        work: Path = root / '.work'
        work.mkdir(exist_ok=True)
        header: str = 'Type,Subtype,Subset,Filter,Genotype,QQ.Field,QQ,Subset.Level,TRUTH.TP\n'
        row: str = 'SNP,*,*,ALL,*,QUAL,*,{},{}\n'
        cases: list[tuple[str, str, str, str, int]] = [
            ('1', '1.0', '9007199254740993', '9007199254740993', 0),
            ('1', '1', '9007199254740993', '9007199254740992', 1),
            ('1', '1', '1', '1.0001', 1),
            ('1', '1.0001', '1', '1', 1),
            ('.', 'NaN', '1', '1', 0),
            ('.', '1', '1', '1', 1),
        ]
        with tempfile.TemporaryDirectory(prefix='extended-contract-', dir=work) as directory:
            left: Path = Path(directory) / 'left.csv'
            right: Path = Path(directory) / 'right.csv'
            def compare() -> int:
                return subprocess.run(
                    [sys.executable, str(root / 'src/sh/compare_extended.py'), str(left), str(right)],
                    capture_output=True, text=True, check=False,
                ).returncode
            for a, b, c, d, expected in cases:
                with self.subTest(values=(a, b, c, d)):
                    left.write_text(header + row.format(a, c))
                    right.write_text(header + row.format(b, d))
                    self.assertEqual(compare(), expected)
            left.write_text(header + row.format('1', '1'))
            right.write_text(header.rstrip() + ',FP.al\n' + row.format('1', '1').rstrip() + ',0\n')
            self.assertNotEqual(compare(), 0)
            right.write_text(header + row.format('1', '1').replace('SNP', 'INDEL'))
            self.assertNotEqual(compare(), 0)
            left.write_text(header + row.format('1', '1').replace('SNP', 'INDEL'))
            right.write_text(header + row.format('1', '1'))
            diagnostic = subprocess.run(
                [sys.executable, str(root / 'src/sh/compare_extended.py'), str(left), str(right)],
                capture_output=True, text=True, check=False,
            )
            self.assertNotEqual(diagnostic.returncode, 0)
            self.assertIn("In 1 but not 2: {('INDEL'", diagnostic.stderr)
            self.assertIn("In 2 but not 1: {('SNP'", diagnostic.stderr)
            left.write_text(header + row.format('1', '1'))
            right.write_text(header + row.format('1', '1') * 2)
            self.assertNotEqual(compare(), 0)
            ratio_header: str = header.rstrip() + ',METRIC.Recall\n'
            left.write_text(ratio_header + row.format('1', '1').rstrip() + ',0.12340\n')
            right.write_text(ratio_header + row.format('1', '1').rstrip() + ',1.234e-1\n')
            self.assertEqual(compare(), 0)
            right.write_text(ratio_header + row.format('1', '1').rstrip() + ',0.1244\n')
            self.assertNotEqual(compare(), 0)


if __name__ == '__main__':
    unittest.main()
