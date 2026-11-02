#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Full-record equality, sorted streaming and format roundtrips."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from Tools.deduplicate import deduplicate, unique_records

HEADER = ('##fileformat=VCFv4.2\n##contig=<ID=c,length=100>\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n')
ROW = 'c\t2\t.\tA\tC\t30\tPASS\t.\tGT\t0/1\n'


class DeduplicateTest(unittest.TestCase):
    def test_full_fields_and_order(self) -> None:
        variants = [ROW, ROW.replace('0/1', '1/1'), ROW.replace('\t30\t', '\t31\t'),
                    ROW.replace('\tC\t', '\tG\t'), ROW.replace('\t.\tA', '\tid\tA'),
                    ROW.replace('PASS', '.'), ROW.replace('\t.\tGT', '\tX=1\tGT'),
                    ROW.replace('GT\t0/1', 'GT:AD\t0/1:1,2')]
        self.assertEqual(list(unique_records(variants + [ROW] + variants)), variants)
        multi = ROW.rstrip() + '\t1/1\n'
        other_sample = ROW.rstrip() + '\t0/1\n'
        self.assertEqual(list(unique_records([multi, other_sample, multi])), [multi, other_sample])
        self.assertEqual(list(unique_records([ROW, ROW.replace('c\t2', 'z\t1')])),
                         [ROW, ROW.replace('c\t2', 'z\t1')])
        for rows in ([ROW, ROW.replace('c\t2', 'c\t1')],
                     [ROW, ROW.replace('c\t2', 'z\t1'), ROW]):
            with self.assertRaises(ValueError):
                list(unique_records(rows))

    def test_formats_and_index_failure(self) -> None:
        with tempfile.TemporaryDirectory(prefix='dedup é ') as directory:
            root = Path(directory)
            plain = root / 'input.vcf'
            plain.write_text(HEADER + ROW + ROW.replace('0/1', '1/1') + ROW)
            for suffix, mode in (('.vcf', 'w'), ('.vcf.gz', 'wz'), ('.bcf', 'wb')):
                path = root / ('out' + suffix)
                if suffix == '.vcf':
                    path.write_bytes(plain.read_bytes())
                else:
                    with pysam.VariantFile(plain) as source:
                        with pysam.VariantFile(path, mode, header=source.header) as dest:
                            for record in source:
                                dest.write(record)
                self.assertEqual(deduplicate(str(path)), (3, 1))
                with pysam.VariantFile(path) as source:
                    self.assertEqual(len(list(source)), 2)
                self.assertEqual(deduplicate(str(path)), (2, 0))
                if suffix != '.vcf':
                    self.assertTrue(Path(str(path) + ('.csi' if suffix == '.bcf' else '.tbi')).exists())
            empty = root / 'empty.vcf'
            empty.write_text(HEADER)
            self.assertEqual(deduplicate(str(empty)), (0, 0))
            broken = root / 'fail.vcf.gz'
            with pysam.BGZFile(str(broken), 'w') as output:
                output.write((HEADER + ROW + ROW).encode())
            original = broken.read_bytes()
            with patch('Tools.deduplicate.runBcftools', side_effect=RuntimeError('index failed')):
                with self.assertRaises(RuntimeError):
                    deduplicate(str(broken))
            self.assertEqual(broken.read_bytes(), original)
            self.assertEqual(list(root.glob('.deduplicate-*')), [])
            plain.write_text(HEADER + ROW + ROW.replace('c\t2', 'c\t1'))
            original_text = plain.read_text()
            with self.assertRaises(ValueError):
                deduplicate(str(plain))
            self.assertEqual(plain.read_text(), original_text)
            self.assertEqual(list(root.glob('.deduplicate-*')), [])


if __name__ == '__main__':
    unittest.main()
