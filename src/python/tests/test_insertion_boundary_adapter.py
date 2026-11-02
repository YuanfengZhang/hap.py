#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify the fixed boundary adapter rejects extra categories, missing records, and fixture changes."""
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


class InsertionBoundaryAdapterTest(unittest.TestCase):
    def test_fixed_contract_rejections(self) -> None:
        root: Path = Path(__file__).resolve().parents[3]
        with tempfile.TemporaryDirectory(dir=root / '.work', prefix='boundary-adapter-') as directory:
            sandbox: Path = Path(directory)
            script: Path = sandbox / 'src/sh/adapt_insertion_boundary_summary.py'
            script.parent.mkdir(parents=True)
            shutil.copyfile(root / 'src/sh/adapt_insertion_boundary_summary.py', script)
            fixture: Path = sandbox / 'example/happy'
            fixture.mkdir(parents=True)
            bed: Path = fixture / 'PG_Conf_hg38-chr21.bed.gz'
            with gzip.open(bed, 'wt') as stream:
                stream.write('chrT\t10\t20\n')
            golden: Path = fixture / 'expected-stratified.summary.csv'
            golden.write_text('Type,Filter,QUERY.FP,QUERY.UNK,QUERY.TOTAL,TRUTH.TP,TRUTH.FN,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score\nINDEL,ALL,1,0,2,1,0,0.5,0,0.6666666666666666\n')
            header: str = '##fileformat=VCFv4.2\n##contig=<ID=chrT,length=100>\n'
            header += '##INFO=<ID=Regions,Number=.,Type=String,Description="Regions">\n'
            for field in ('GT', 'BD', 'BK', 'BVT', 'BLT'):
                header += f'##FORMAT=<ID={field},Number=1,Type=String,Description="{field}">\n'
            header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n'
            record: str = 'chrT\t10\t.\tC\tCT\t40\tPASS\t.\tGT:BD:BK:BVT:BLT\t./.:.:.:.:.\t0/1:UNK:lm:INDEL:het\n'
            vcf: Path = sandbox / 'actual.vcf'
            vcf.write_text(header + record)
            rows: list = [[['chrT', 10, 'C', ['CT']], ['PASS'], {
                'TRUTH': {'GT': [None, None], 'BD': '.', 'BK': '.', 'BVT': '.', 'BLT': '.'},
                'QUERY': {'GT': [0, 1], 'BD': 'UNK', 'BK': 'lm', 'BVT': 'INDEL', 'BLT': 'het'},
            }]]
            contract: dict = {
                'version': 'synthetic-known-answer', 'insertion_keys': [['chrT', 10, 'C', ['CT']]],
                'record_count': 1,
                'expected_records': 'insertion-boundary-expected-v1.jsonl.gz',
                'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (bed, golden)},
            }
            contract_path: Path = sandbox / 'src/data/insertion-boundary-summary-v1.json'
            contract_path.parent.mkdir()
            contract_path.write_text(json.dumps(contract))
            with gzip.open(contract_path.parent / contract['expected_records'], 'wt') as stream:
                for row in rows:
                    stream.write(json.dumps(row) + '\n')

            def invoke(name: str) -> subprocess.CompletedProcess:
                return subprocess.run([sys.executable, '-I', '-O', str(script), '--actual-vcf', str(vcf),
                                       '--output', str(sandbox / (name + '.summary.csv'))],
                                      capture_output=True, text=True, check=False)

            positive = invoke('positive')
            self.assertEqual(positive.returncode, 0, positive.stderr)
            for name, contents in [('extra-classification', header + record.replace('UNK', 'FP')),
                                    ('missing-record', header),
                                    ('extra-record', header + record + record.replace('chrT\t10', 'chrT\t11')),
                                    ('filter-changed', header + record.replace('PASS', '.')),
                                    ('gt-changed', header + record.replace('0/1:UNK', '1/1:UNK'))]:
                vcf.write_text(contents)
                rejected = invoke(name)
                self.assertNotEqual(rejected.returncode, 0)
                self.assertIn('Unexpected record or classification change', rejected.stderr)
            vcf.write_text(header + record)
            golden.write_text(golden.read_text() + '\n')
            rejected = invoke('fixture-changed')
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn('Fixed input changed', rejected.stderr)
            contract['input_sha256'][golden.name] = hashlib.sha256(golden.read_bytes()).hexdigest()
            contract_path.write_text(json.dumps(contract))
            with gzip.open(bed, 'wt') as stream:
                stream.write('chrT\t9\t20\n')
            rejected = invoke('bed-changed')
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn('Fixed input changed', rejected.stderr)


if __name__ == '__main__':
    unittest.main()
