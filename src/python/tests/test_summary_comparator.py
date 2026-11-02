#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Positive and negative regressions for explicit schema allowlists and exact scientific counts."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class SummaryComparatorTest(unittest.TestCase):
    def test_versioned_schema(self) -> None:
        root: Path = Path(__file__).resolve().parents[3]
        (root / '.work').mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=root / '.work', prefix='summary-contract-') as directory:
            work: Path = Path(directory)
            actual: Path = work / 'case.summary.csv'
            expected: Path = work / 'golden.csv'
            extended: Path = work / 'case.extended.csv'
            old: dict[str, str] = {'Type': 'SNP', 'Filter': 'ALL', 'TRUTH.TP': '1',
                                   'TRUTH.FN': '1', 'QUERY.TOTAL': '2', 'QUERY.FP': '1'}
            new: dict[str, str] = dict(old, **{'FP.al': '1', 'METRIC.F1_Score': '0.5'})
            ext: dict[str, str] = dict(new, **{'QUERY.TP': '1', 'Subtype': '*',
                                              'Subset': '*', 'Genotype': '*', 'QQ': '*', 'QQ.Field': 'QUAL'})

            def write(path: Path, row: dict[str, str], duplicate: bool = False) -> None:
                with path.open('w', newline='') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(row))
                    writer.writeheader()
                    writer.writerow(row)
                    if duplicate:
                        writer.writerow(row)

            def compare(versioned: bool = True) -> int:
                command: list[str] = [sys.executable, '-I', str(root / 'src/sh/compare_summaries.py'),
                                      str(actual), str(expected)]
                if versioned:
                    command += ['--schema', 'ga4gh-summary-v2']
                return subprocess.run(command, capture_output=True, check=False).returncode

            write(expected, old)
            write(actual, new)
            write(extended, ext)
            self.assertEqual(compare(), 0)
            self.assertNotEqual(compare(False), 0)
            for field, value in (('unknown', '0'), ('FP.al', '0.5'),
                                 ('METRIC.F1_Score', '0.6'), ('QUERY.FP', '2')):
                with self.subTest(field=field):
                    write(actual, dict(new, **{field: value}))
                    self.assertNotEqual(compare(), 0)
            write(actual, new)
            missing_column: dict[str, str] = dict(new)
            del missing_column['QUERY.FP']
            write(actual, missing_column)
            self.assertNotEqual(compare(), 0)
            write(actual, new)
            extended.unlink()
            self.assertNotEqual(compare(), 0)
            write(extended, dict(ext, **{'QQ.Field': 'GQ'}))
            self.assertNotEqual(compare(), 0)
            write(extended, ext, duplicate=True)
            self.assertNotEqual(compare(), 0)
            write(extended, ext)
            write(actual, dict(new, **{'TRUTH.TP': '9007199254740993'}))
            write(expected, dict(old, **{'TRUTH.TP': '9007199254740992'}))
            self.assertNotEqual(compare(False), 0)
            for missing in ('', '.', 'NaN'):
                for swapped in (False, True):
                    with self.subTest(missing=missing, swapped=swapped):
                        a, b = (missing, '1') if not swapped else ('1', missing)
                        write(actual, dict(old, **{'TRUTH.TP': a}))
                        write(expected, dict(old, **{'TRUTH.TP': b}))
                        self.assertNotEqual(compare(False), 0)
            exact: dict[str, str] = dict(old, **{'TRUTH.TP': '9007199254740993'})
            write(actual, exact)
            write(expected, dict(exact, **{'TRUTH.TP': '9007199254740992'}))
            self.assertNotEqual(compare(False), 0)


if __name__ == '__main__':
    unittest.main()
