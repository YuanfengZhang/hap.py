#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Characterize quarantine records, reasons and strict parse failures."""
import gzip
from pathlib import Path
import sys
import tempfile
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Tools.invalid import check_outputs, record_errors, split_invalid

HEADER = ('##fileformat=VCFv4.2\n##contig=<ID=c,length=100>\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
          '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n')
GOOD = 'c\t1\t.\tA\tC\t30\tPASS\t.\tGT:AD\t0|1:8,2\n'
BAD = 'c\t2\t.\tA\tC\t30\tPASS\t.\tGT:AD\t0/2:1,2\n'


class InvalidTest(unittest.TestCase):
    def test_rules(self) -> None:
        for gt in ('0/1', '1', '.', './.', '.|1'):
            self.assertEqual(record_errors(GOOD.replace('0|1', gt), ['S']), [])
        self.assertEqual(record_errors(BAD, ['S'])[0][1], 'GT_ALLELE')
        self.assertEqual(record_errors(GOOD.replace('0|1', '0/1/1'), ['S'])[0][1], 'GT_PLOIDY')
        self.assertEqual(record_errors(GOOD.replace('8,2', '8,2,3'), ['S'])[0][1], 'AD_CAPACITY')
        self.assertEqual(record_errors(GOOD.replace('\tC\t', '\tC,C\t').replace('0|1:8,2', '2/1:0,8,13'), ['S']), [])
        multi = GOOD.rstrip() + '\t0/2:1,2\n'
        self.assertEqual(record_errors(multi, ['S', 'Q'])[0][0], 'Q')

    def test_formats_and_raw_preservation(self) -> None:
        with tempfile.TemporaryDirectory(prefix='invalid é space ') as directory:
            root = Path(directory)
            for kind in ('vcf', 'vcf.gz', 'bcf'):
                source = root / ('input.' + kind)
                if kind == 'vcf':
                    source.write_text(HEADER + GOOD + BAD)
                elif kind == 'vcf.gz':
                    with gzip.open(source, 'wt') as stream:
                        stream.write(HEADER + GOOD + BAD)
                else:
                    with pysam.VariantFile(root / 'input.vcf') as src:
                        with pysam.VariantFile(source, 'wb', header=src.header) as dest:
                            for record in src:
                                dest.write(record)
                output, valid = str(root / (kind + '.invalid.gz')), str(root / (kind + '.valid.gz'))
                result = split_invalid(str(source), output, valid)
                self.assertEqual((result.total, result.invalid), (2, 1))
                with gzip.open(output, 'rt') as stream:
                    rejected = stream.read()
                self.assertTrue(rejected.endswith(BAD))
                if kind != 'bcf':
                    self.assertEqual(rejected, HEADER + BAD)
                with gzip.open(valid, 'rt') as stream:
                    self.assertTrue(stream.read().endswith(GOOD))
                self.assertIn('2\tc\t2\tS\tGT_ALLELE', Path(output + '.reasons.tsv').read_text())
                with self.assertRaises(FileExistsError):
                    split_invalid(str(source), output, valid)

    def test_header_and_stream_errors_stay_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cases = {
                'contig': HEADER + GOOD.replace('c\t1', 'absent\t1'),
                'format': HEADER + GOOD.replace('GT:AD', 'GT:UNDECLARED'),
                'syntax': HEADER + GOOD.replace('0|1', '0/q'),
                'extra-sample': HEADER + GOOD.rstrip() + '\t0/1\n',
            }
            for label, content in cases.items():
                source = root / (label + '.vcf')
                output = root / (label + '.invalid.gz')
                source.write_text(content)
                with self.assertRaises((ValueError, OSError)):
                    split_invalid(str(source), str(output), str(root / 'valid.gz'))
                self.assertFalse(output.exists())
            source = root / 'truncated.vcf.gz'
            source.write_bytes(gzip.compress((HEADER + GOOD).encode())[:-6])
            with self.assertRaises((ValueError, OSError, EOFError)):
                split_invalid(str(source), str(root / 'broken.invalid.gz'), str(root / 'valid.gz'))
            self.assertFalse((root / 'broken.invalid.gz').exists())
            self.assertEqual(list(root.glob('.invalid-pending-*')), [])
            self.assertEqual(list(root.glob('invalid-reader-*')), [])

    def test_empty_clean_and_failure(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for label, rows in [('good', GOOD), ('bad', BAD), ('malformed', GOOD + 'c\tbad\n')]:
                source, output, valid = [str(root / (label + suffix)) for suffix in ('.vcf', '.invalid.gz', '.valid.gz')]
                Path(source).write_text(HEADER + rows)
                if label == 'good':
                    self.assertEqual(split_invalid(source, output, valid).invalid, 0)
                    with gzip.open(output, 'rt') as stream:
                        self.assertEqual(stream.read(), HEADER)
                else:
                    with self.assertRaises((ValueError, OSError)):
                        split_invalid(source, output, valid)
                    self.assertEqual(Path(output).exists(), label == 'bad')
            self.assertEqual(list(root.glob('.invalid-pending-*')), [])
            with self.assertRaises(ValueError):
                check_outputs([str(root / 'good.vcf')], [str(root / 'good.vcf')])
            dangling = root / 'dangling.gz'
            dangling.symlink_to(root / 'absent.gz')
            with self.assertRaises(FileExistsError):
                check_outputs([str(dangling)], [])


if __name__ == '__main__':
    unittest.main()
