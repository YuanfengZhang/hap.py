#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Keep caller registries, field ordering and tag-specific dispatch compatible."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import Somatic


class FeatureRegistryTest(unittest.TestCase):
    def test_registered_callers_and_truth_fields(self) -> None:
        expected = {
            'hcc.mutect.snv': 'extractMutectSNVFeatures',
            'hcc.mutect.indel': 'extractMutectIndelFeatures',
            'hcc.varscan2.snv': 'extractVarscan2SNVFeatures',
            'hcc.varscan2.indel': 'extractVarscan2IndelFeatures',
            'hcc.pisces.snv': 'extractPiscesSNVFeatures',
            'hcc.pisces.indel': 'extractPiscesIndelFeatures',
        }
        for name, extractor in expected.items():
            with self.subTest(name=name):
                feature = Somatic.FeatureSet.make(name)
                feature.setChrDepths({'chr1': 12.})
                with patch.object(Somatic, extractor, return_value='caller') as collect:
                    self.assertEqual(feature.collect('query.vcf', 'FP'), 'caller')
                    collect.assert_called_once_with('query.vcf', 'FP', {'chr1': 12.})
                with patch.object(Somatic.GenericFeatures, 'collectFeatures', return_value='truth') as collect:
                    self.assertEqual(feature.collect('truth.vcf', 'TP'), 'truth')
                    fields = collect.call_args.args[2]
                    self.assertEqual(fields[:5], ['CHROM', 'POS', 'REF', 'ALT', 'QUAL'])
                    self.assertEqual(fields[-5:], ['I.T_ALT_RATE', 'I.DP_normal', 'I.DP_tumor', 'I.tag', 'I.count'])

    def test_list_conversion(self) -> None:
        self.assertEqual(Somatic.GenericFeatures.processValue(('ALT', ['A', 'C'])), 'A,C')
        self.assertIsNone(Somatic.GenericFeatures.processValue(('QUAL', None)))


if __name__ == '__main__':
    unittest.main()
