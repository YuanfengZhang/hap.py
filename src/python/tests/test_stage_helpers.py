#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exercise report-region resolution and somatic depth aggregation seams."""
from argparse import Namespace
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd
import qfy
import som


class StageHelpersTest(unittest.TestCase):
    def test_region_resolution_and_duplicate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / 'regions.bed'
            bed.write_text('chr1\t0\t10\n')
            strat = root / 'strata.tsv'
            strat.write_text('relative\tregions.bed\n\n')
            args = Namespace(fp_bedfile=str(bed), strat_tsv=str(strat),
                             strat_regions=[f'explicit:{bed}'])
            result = qfy.resolve_regions(args)
            self.assertEqual(list(result), ['CONF', 'relative', 'explicit'])
            self.assertEqual(result['relative'], str(bed))
            strat.write_text('CONF\tregions.bed\n')
            with self.assertRaisesRegex(Exception, 'Duplicate stratification'):
                qfy.resolve_regions(args)

    def test_missing_region(self) -> None:
        args = Namespace(fp_bedfile=None, strat_tsv=None,
                         strat_regions=['missing:/nonexistent/happy-regions.bed'])
        with self.assertRaisesRegex(Exception, 'not found'):
            qfy.resolve_regions(args)

    def test_coverage_aggregation(self) -> None:
        self.assertIsNone(som.coverage_normalization([]))
        frames = [pd.DataFrame({'CHROM': ['chr1', 'chr2'], 'COVERAGE': [10., 4.]}),
                  pd.DataFrame({'CHROM': ['chr1'], 'COVERAGE': [20.]})]
        with patch.object(som, 'bamStats', side_effect=frames) as collect:
            result = som.coverage_normalization(['first.bam', 'second.bam'])
        self.assertEqual(result, {'chr1': 45., 'chr2': 12.})
        self.assertEqual(collect.call_count, 2)


if __name__ == '__main__':
    unittest.main()
