#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Characterize the hap.py CLI independently of comparison execution."""
import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src/python'))
SPEC = importlib.util.spec_from_file_location('happy_cli', ROOT / 'src/python/hap.py')
CLI = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLI)


class ParserTest(unittest.TestCase):
    def test_defaults(self) -> None:
        with patch.object(CLI.Tools, 'has_sge', False):
            args = CLI.build_parser().parse_args(['truth.vcf', 'query.vcf'])
        self.assertEqual(args._vcfs, ['truth.vcf', 'query.vcf'])
        self.assertEqual(args.engine, 'xcmp')
        self.assertFalse(args.output_invalid)
        self.assertFalse(args.keep_duplicates)
        self.assertTrue(args.delete_scratch)
        self.assertIsNone(args.scratch_prefix)
        self.assertIsNone(args.reports_prefix)
        self.assertEqual(args.engine_scmp_distance, 30)

    def test_explicit_options(self) -> None:
        args = CLI.build_parser().parse_args([
            'truth.vcf', 'query.vcf', '--engine', 'scmp-distance',
            '--keep-scratch', '--scratch-prefix', 'scratch space',
            '--threads', '2', '--no-roc', '-o', 'report', '--output-invalid', '--keep-duplicates',
        ])
        self.assertEqual(args.engine, 'scmp-distance')
        self.assertFalse(args.delete_scratch)
        self.assertFalse(args.do_roc)
        self.assertEqual(args.roc, 'QUAL')
        self.assertEqual(args.threads, 2)
        self.assertTrue(args.output_invalid)
        self.assertTrue(args.keep_duplicates)
        self.assertEqual(args.scratch_prefix, 'scratch space')
        self.assertEqual(args.reports_prefix, 'report')

    def test_build_has_no_execution_side_effects(self) -> None:
        with patch.object(CLI.Tools, 'init') as initialize, \
                patch.object(CLI.qfy, 'quantify') as quantify:
            CLI.build_parser()
        initialize.assert_not_called()
        quantify.assert_not_called()


if __name__ == '__main__':
    unittest.main()
