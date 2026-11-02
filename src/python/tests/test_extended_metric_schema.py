#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Allow zero-denominator extended metrics only under an explicit versioned contract."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class ExtendedMetricSchemaTest(unittest.TestCase):
    def test_zero_denominator_contract(self) -> None:
        root: Path = Path(__file__).resolve().parents[3]
        base: dict[str, str] = dict(zip(
            ('Type', 'Subtype', 'Subset', 'Filter', 'Genotype', 'QQ.Field', 'QQ'),
            ('INDEL', '*', '*', 'ALL', '*', 'QUAL', '*')))
        base.update({'TRUTH.TP': '0', 'TRUTH.FN': '1', 'QUERY.TP': '0', 'QUERY.FP': '1',
                     'QUERY.TOTAL': '1', 'QUERY.UNK': '0', 'METRIC.Recall': '0',
                     'METRIC.Precision': '0', 'METRIC.F1_Score': '0', 'METRIC.Frac_NA': '0'})
        with tempfile.TemporaryDirectory(dir=root / '.work', prefix='extended-metrics-') as directory:
            actual: Path = Path(directory) / 'actual.csv'
            expected: Path = Path(directory) / 'expected.csv'

            def write(path: Path, row: dict[str, str]) -> None:
                with path.open('w') as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(row))
                    writer.writeheader()
                    writer.writerow(row)

            def compare(versioned: bool, optimized: bool = False) -> int:
                command: list[str] = [sys.executable, '-I', str(root / 'src/sh/compare_extended.py'),
                                      str(actual), str(expected)]
                if versioned:
                    command += ['--schema', 'ga4gh-metrics-v2']
                if optimized:
                    command.insert(1, '-O')
                return subprocess.run(command, capture_output=True, check=False).returncode

            write(actual, base)
            write(expected, dict(base, **{'METRIC.F1_Score': ''}))
            self.assertNotEqual(compare(False), 0)
            self.assertEqual(compare(True), 0)
            write(actual, dict(base, **{'METRIC.F1_Score': '1'}))
            self.assertNotEqual(compare(True), 0)
            write(actual, dict(base, **{'QUERY.FP': '2'}))
            self.assertNotEqual(compare(True), 0)
            positive: dict[str, str] = dict(base, **{'TRUTH.TP': '1', 'TRUTH.FN': '0',
                'QUERY.TP': '1', 'QUERY.FP': '0', 'METRIC.Recall': '1',
                'METRIC.Precision': '1', 'METRIC.F1_Score': '1'})
            write(actual, positive)
            write(expected, dict(positive, **{'METRIC.F1_Score': ''}))
            self.assertNotEqual(compare(True), 0)
            write(expected, positive)
            write(actual, dict(positive, **{'METRIC.Precision': '0.9'}))
            self.assertNotEqual(compare(True), 0)
            write(expected, dict(positive, **{'METRIC.Precision': '0.9'}))
            self.assertNotEqual(compare(True, optimized=True), 0)
            for value in ('NaN', 'Infinity', '-1'):
                invalid: dict[str, str] = dict(base, **{'QUERY.TP': value})
                write(actual, invalid)
                write(expected, invalid)
                self.assertNotEqual(compare(True, optimized=True), 0)
            write(actual, base)
            write(expected, dict(base, **{'UNKNOWN': '0'}))
            self.assertNotEqual(compare(True), 0)
            missing: dict[str, str] = dict(base)
            del missing['QUERY.TP']
            write(actual, missing)
            write(expected, missing)
            self.assertNotEqual(compare(True), 0)
            write(actual, base)
            write(expected, base)
            with actual.open('a') as stream:
                csv.DictWriter(stream, fieldnames=list(base)).writerow(base)
            self.assertNotEqual(compare(True), 0)


if __name__ == '__main__':
    unittest.main()
