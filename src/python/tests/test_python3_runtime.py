#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python 3 runtime and independent benchmark metrics regressions."""
import ast
import gzip
import importlib
import io
import multiprocessing
import os
from pathlib import Path
import pickle
import sys
import subprocess
import shutil
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'src/python'))
import numpy as np
import pandas as pd
import pysam
import Tools
from Tools.bcftools import runShellCommand, concatenateParts, parseStats
from Tools.parallel import runParallel, unpickleSequentially
from Haplo.happyroc import roc, RESULT_ALLCOLUMNS


def square(value: int) -> int:
    return value * value


def fail(value: int) -> int:
    raise ValueError('worker failure')


class RuntimeTests(unittest.TestCase):
    def test_syntax_and_imports(self):
        for path in (ROOT / 'src/python').rglob('*.py'):
            ast.parse(path.read_text(), filename=str(path))
        for path in (ROOT / 'src/sh').glob('*.py'):
            ast.parse(path.read_text(), filename=str(path))
        ast.parse((ROOT / '.ycm_extra_conf.py').read_text())
        template = (ROOT / 'src/python/Haplo/version.py.in').read_text()
        ast.parse(template.replace('@HAPLOTYPES_VERSION@', 'test').replace('@SGE_AVAILABLE@', 'False').replace('@VCFEVAL_AVAILABLE@', 'False'))
        for package in ('Tools', 'Haplo', 'Somatic'):
            for path in (ROOT / 'src/python' / package).glob('*.py'):
                importlib.import_module(package + '.' + path.stem)

    def test_pipeline(self):
        py = sys.executable
        self.assertEqual(runShellCommand(py, '-c', 'print("峰 space")'), '峰 space\n')
        self.assertEqual(runShellCommand(py, '-c', 'import sys;sys.stdout.buffer.write(bytes(range(256)))', binary=True), bytes(range(256)))
        with self.assertRaises(UnicodeDecodeError):
            runShellCommand(py, '-c', 'import sys;sys.stdout.buffer.write(b"\\xff")')
        for first, last in [(3, 0), (0, 4)]:
            with self.assertRaises(RuntimeError):
                runShellCommand(py, '-c', f'import sys;sys.stderr.write("x"*100000);sys.exit({first})', '|', py, '-c', f'import sys;sys.stdin.buffer.read();sys.exit({last})')
        self.assertEqual(runShellCommand(py, '-c', 'print("abc")', '|', py, '-c', 'import sys; print(sys.stdin.read().upper(),end="")'), 'ABC\n')
        with tempfile.TemporaryDirectory(prefix='峰 space ') as directory:
            path = Path(directory) / '峰 file'
            path.write_text('hello')
            self.assertEqual(runShellCommand(py, '-c', 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).read_text())', str(path)), 'hello\n')

    def test_reference_and_bgzf(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory)/'a.fa', Path(directory)/'b.fa'
            first.write_text('>custom\nACGT\n')
            second.write_text('>other\nAAA\n')
            pysam.faidx(str(first))
            pysam.faidx(str(second))
            with patch.dict(os.environ, {'HGREF': str(first), 'HG19': str(second)}):
                self.assertEqual(Tools.defaultReference(), str(first))
                self.assertEqual(Tools.defaultReference(str(second)), str(second))
                with self.assertRaises(ValueError):
                    Tools.defaultReference(str(first)+'missing')
                header = io.StringIO()
                Tools.writeVCFHeader(header)
                self.assertIn('##contig=<ID=custom,length=4>', header.getvalue())
                self.assertNotIn('hg19', header.getvalue())
            target = str(Path(directory)/'峰.vcf.gz')
            with Tools.BGZipFile(target) as stream:
                stream.write('峰\n')
                stream.write(b'bytes\n')
            with gzip.open(target, 'rt') as stream:
                self.assertEqual(stream.read(), '峰\nbytes\n')

    def test_pools_and_pickle(self):
        for method in ('spawn', 'forkserver'):
            if method not in multiprocessing.get_all_start_methods():
                continue
            with multiprocessing.get_context(method).Pool(2) as pool:
                self.assertEqual(runParallel(pool, square, range(5)), [0, 1, 4, 9, 16])
                with self.assertRaisesRegex(ValueError, 'worker failure'):
                    runParallel(pool, fail, [1])
        self.assertEqual(runParallel(None, square, range(5)), [0, 1, 4, 9, 16])
        with self.assertRaises(ValueError):
            runParallel(None, fail, [1])
        with tempfile.NamedTemporaryFile(delete=False) as stream:
            pickle.dump(['峰', 2], stream)
            name = stream.name
        self.assertEqual(list(unpickleSequentially([name])), ['峰', 2])
        self.assertFalse(Path(name).exists())

    def test_concatenate_ten_parts(self):
        calls = []
        def fake(*args):
            calls.append(args)
            if args[0] == 'concat':
                output = args[args.index('-o')+1]
                inputs = args[args.index('-o')+2:]
                Path(output).write_bytes(b''.join(Path(p).read_bytes() for p in inputs))
            return ''
        with tempfile.TemporaryDirectory() as directory:
            parts = [str(Path(directory)/f'{i}.bcf') for i in range(12)]
            for i, part in enumerate(parts):
                Path(part).write_bytes(bytes([i]))
            output = str(Path(directory)/'out.bcf')
            with patch('Tools.bcftools.runBcftools', side_effect=fake):
                concatenateParts(output, *parts)
            self.assertEqual(Path(output).read_bytes(), bytes(range(12)))
            self.assertGreater(len([c for c in calls if c[0]=='concat']), 1)
        self.assertEqual(parseStats('SN\t0\tnumber of records:\t12\n').iloc[0]['count'], 12)

    def test_roc_strata_and_zero_denominators(self):
        rows=[]
        for subset, counts in [('*',(2,1,1,1,0)), ('other',(1,0,3,1,0)), ('unknown',(0,0,0,0,5)), ('empty',(0,0,0,0,0)), ('zero',(0,2,0,3,0))]:
            tt, fn, qt, fp, unk = counts
            row = dict.fromkeys(RESULT_ALLCOLUMNS, 0)
            row.update({'Type':'SNP','Subtype':'*','Subset':subset,'Filter':'ALL','Genotype':'*','QQ.Field':'QUAL','QQ':'10', 'TRUTH.TP':tt,'TRUTH.FN':fn,'TRUTH.TOTAL':tt+fn,'QUERY.TP':qt,'QUERY.FP':fp,'QUERY.UNK':unk,'QUERY.TOTAL':qt+fp+unk})
            rows.append(row)
        with tempfile.TemporaryDirectory() as directory:
            input_path=Path(directory)/'roc.tsv'
            pd.DataFrame(rows).to_csv(input_path, sep='\t', index=False)
            result=roc(input_path, str(Path(directory)/'result'), total_region_size=400)['all'].set_index('Subset')
            self.assertEqual(result.loc['*','Subset.Size'],400)
            self.assertAlmostEqual(result.loc['*','METRIC.Recall'], 2/3)
            self.assertEqual(result.loc['*','METRIC.Precision'], .5)
            self.assertAlmostEqual(result.loc['*','METRIC.F1_Score'], 4/7)
            self.assertEqual(result.loc['other','QUERY.TP'], 3)
            for subset in ('unknown','empty','zero'):
                self.assertEqual(result.loc[subset,'METRIC.Precision'], 0)
                self.assertEqual(result.loc[subset,'METRIC.F1_Score'], 0)
            reread=pd.read_csv(Path(directory)/'result.all.csv.gz').set_index('Subset')
            np.testing.assert_allclose(reread['METRIC.Precision'], result['METRIC.Precision'])

    def test_subset_integer_metadata(self):
        from Haplo.happyroc import _postprocessRocData
        rows=[]
        for value in ('NA', '.', '', '402.000000', '9007199254740993', '-9223372036854775808', '9223372036854775807'):
            row=dict.fromkeys(RESULT_ALLCOLUMNS,0)
            row.update({'Type':'SNP','Subset':str(value),'Subtype':'*','Filter':'ALL','Genotype':'*','QQ.Field':'QUAL','QQ':'*','Subset.Size':value})
            rows.append(row)
        result=_postprocessRocData(pd.DataFrame(rows)).set_index('Subset')
        self.assertTrue(pd.isna(result.loc['NA','Subset.Size']))
        self.assertEqual(result.loc['402.000000','Subset.Size'],402)
        self.assertEqual(result.loc['-9223372036854775808','Subset.Size'],-9223372036854775808)
        self.assertEqual(result.loc['9007199254740993','Subset.Size'],9007199254740993)
        self.assertEqual(result.loc['9223372036854775807','Subset.Size'],9223372036854775807)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'subset.csv'
            result.to_csv(path)
            restored=pd.read_csv(path, dtype={'Subset.Size':'Int64'},keep_default_na=False,na_values=['']).set_index('Subset')
            self.assertEqual(restored.loc['9007199254740993','Subset.Size'],9007199254740993)
        from Tools.metric import dataframeToMetricsTable
        import json
        payload = json.loads(json.dumps(dataframeToMetricsTable('subset', result), allow_nan=False))
        column = next(item for item in payload['data'] if item['id']=='Subset.Size')
        self.assertEqual(column['type'], 'int64')
        self.assertIn(9007199254740993, column['values'])
        self.assertIn(None, column['values'])
        self.assertTrue(all(value is None or type(value) is int for value in column['values']))
        thresholds = [dict(rows[3], **{'QQ': value}) for value in ('40', '2', '*', '10')]
        ordered = _postprocessRocData(pd.DataFrame(thresholds))
        self.assertEqual(ordered['QQ'].tolist(), ['*', '2', '10', '40'])
        for bad in ('invalid', '1.5', '9223372036854775808', '-9223372036854775809', 'Infinity', '-Infinity', 'NaN'):
            invalid=pd.DataFrame([dict(rows[0], **{'Subset.Size':bad})])
            with self.assertRaises((ValueError,OverflowError,TypeError)):
                _postprocessRocData(invalid)

    def test_bed_json_and_session_info(self):
        import json
        from Tools.bedintervaltree import BedIntervalTree
        from Tools.bcftools import bedOverlapCheck
        from Tools.metric import dataframeToMetricsTable
        from Tools.sessioninfo import sessionInfo
        table = dataframeToMetricsTable('test', pd.DataFrame({'n':[1,2], 'f':[np.nan,1.5]}))
        json.dumps(table, allow_nan=False)
        info = sessionInfo()
        self.assertTrue(info['pip_packages'])
        json.dumps(info)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'regions.bed.gz'
            with gzip.open(path,'wt') as stream:
                stream.write('# regions\ncustom\t0\t10\ncustom\t10\t20\n')
            self.assertEqual(bedOverlapCheck(str(path)),0)
            tree=BedIntervalTree()
            tree.addFromBed(str(path),'CONF')
            self.assertEqual(tree.countbases(label='CONF'),20)
            self.assertEqual(len(tree.intersect('custom',1,2)),1)
            with gzip.open(path,'wt') as stream:
                stream.write('custom\t0\t10\ncustom\t9\t20\n')
            self.assertEqual(bedOverlapCheck(str(path)),1)

    def test_bcftools_real_roundtrip(self):
        from Tools.bcftools import preprocessVCF, runBcftools
        binary = Path(shutil.which('bcftools') or (Path(sys.executable).parent / 'bcftools'))
        if not binary.exists():
            self.skipTest('bcftools is not installed alongside this interpreter')
        with patch.dict(os.environ, {'PATH': str(binary.parent)+os.pathsep+os.environ.get('PATH','')}):
            with tempfile.TemporaryDirectory(prefix='峰 space ') as directory:
                directory=Path(directory)
                reference=directory/'reference.fa'
                reference.write_text('>custom\n'+'A'*100+'\n')
                pysam.faidx(str(reference))
                header='##fileformat=VCFv4.2\n##contig=<ID=custom,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
                source=directory/'in.vcf'
                source.write_text(header+'custom\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n')
                output=directory/'out.vcf.gz'
                preprocessVCF(str(source),str(output),chrprefix=False,norm=True,reference=str(reference),num_threads=1)
                self.assertTrue(Path(str(output)+'.tbi').exists())
                bcf=directory/'out.bcf'
                runBcftools('view','-Ob','-o',str(bcf),str(output))
                runBcftools('index',str(bcf))
                self.assertTrue(Path(str(bcf)+'.csi').exists())
                self.assertEqual(runBcftools('view','-H',str(bcf)),runBcftools('view','-H',str(output)))
                self.assertTrue(runBcftools('view','-Ob',str(bcf),binary=True).startswith(b'\x1f\x8b'))
                parts=[]
                for i in range(12):
                    part=directory/f'part{i}.vcf'
                    part.write_text(header+f'custom\t{i+1}\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n')
                    compressed=str(part)+'.gz'
                    runBcftools('view','-Oz','-o',compressed,str(part))
                    parts.append(compressed)
                combined=str(directory/'combined.vcf.gz')
                concatenateParts(combined,*parts)
                self.assertEqual(len(runBcftools('view','-H',combined).splitlines()),12)
                source.write_text(header+'custom\t3\t.\tA\tC,<NON_REF>\t10\tPASS\t.\tGT\t0/1\n')
                converted=str(directory/'converted.vcf.gz')
                preprocessVCF(str(source),converted,chrprefix=False,reference=str(reference),convert_gvcf=True,num_threads=1)
                self.assertNotIn('<NON_REF>',runBcftools('view','-H',converted))

    def test_cli_help(self):
        for name in ('hap.py', 'pre.py', 'qfy.py', 'som.py', 'ftx.py', 'cnx.py', 'bamstats.py', 'ovc.py'):
            with self.subTest(cli=name):
                result = subprocess.run([sys.executable, '-I', '-B', str(ROOT/'src/python'/name), '--help'], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout.lower())

    def test_somatic_summary_extended(self):
        import som
        args = SimpleNamespace(features='generic.snv', af_strat_binsize=[.5], af_strat_truth='TAF', af_strat_query='QAF')
        rows = [
            {'REF.truth':'A', 'REF':'A', 'ALT':'C', 'tag':'TP', 'FILTER':'PASS', 'TAF':.2, 'QAF':.7},
            {'REF.truth':'A', 'REF':'A', 'ALT':'C', 'tag':'TP', 'FILTER':'LowQual', 'TAF':.2, 'QAF':.7},
            {'REF.truth':'A', 'REF':None, 'ALT':None, 'tag':'FN', 'FILTER':None, 'TAF':.2, 'QAF':np.nan},
            {'REF.truth':None, 'REF':'A', 'ALT':'C', 'tag':'FP', 'FILTER':'PASS', 'TAF':np.nan, 'QAF':.7},
            {'REF.truth':None, 'REF':'A', 'ALT':'C', 'tag':'UNK', 'FILTER':'PASS', 'TAF':np.nan, 'QAF':.7},
        ]
        frame = pd.DataFrame(rows)
        summary = som.summary_from_featuretable(frame, args).set_index('Filter')
        self.assertEqual(summary.loc['PASS','TRUTH.TP'], 2)
        self.assertEqual((summary.loc['PASS','QUERY.TOTAL'] - summary.loc['PASS','QUERY.FP'] - summary.loc['PASS','QUERY.UNK']), 1)
        self.assertEqual(summary.loc['PASS','QUERY.TOTAL'], 3)
        self.assertEqual(summary.loc['PASS','METRIC.Precision'], .5)
        self.assertEqual(summary.loc['PASS','METRIC.Recall'], .6667)
        self.assertAlmostEqual(summary.loc['PASS','METRIC.F1_Score'], .5714)
        self.assertEqual((summary.loc['ALL','QUERY.TOTAL'] - summary.loc['ALL','QUERY.FP'] - summary.loc['ALL','QUERY.UNK']), 2)
        extended = som.extended_from_featuretable(frame, args).set_index(['Subset','Filter'])
        self.assertNotIn('QUERY.TP', extended.columns)
        self.assertNotIn('QUERY.TP', summary.columns)
        extended['QUERY.TP'] = extended['QUERY.TOTAL'] - extended['QUERY.FP'] - extended['QUERY.UNK']
        self.assertEqual(extended.loc[('[0.00,0.50)','PASS'),'TRUTH.TP'], 2)
        self.assertEqual(extended.loc[('[0.00,0.50)','PASS'),'QUERY.TP'], 0)
        self.assertEqual(extended.loc[('[0.50,1.00]','PASS'),'TRUTH.TP'], 0)
        self.assertEqual(extended.loc[('[0.50,1.00]','PASS'),'QUERY.TP'], 1)
        for data in (frame.iloc[0:0], frame[frame['tag']=='UNK']):
            for output in (som.summary_from_featuretable(data,args), som.extended_from_featuretable(data,args)):
                self.assertTrue((output['METRIC.Precision']==0).all())
                self.assertTrue((output['METRIC.F1_Score']==0).all())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'summary.csv'
            summary.to_csv(path)
            restored = pd.read_csv(path).set_index('Filter')
            np.testing.assert_array_equal(restored['QUERY.TOTAL'], summary['QUERY.TOTAL'])

    def test_comparators_positive_and_negative(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = Path(directory)/'a.csv', Path(directory)/'b.csv'
            row = {'Type':'SNP','Subtype':'*','Subset':'*','Filter':'ALL','Genotype':'*','QQ.Field':'QUAL','QQ':'10','TRUTH.TP':2,'QUERY.TP':1}
            for script in ('compare_summaries.py','compare_extended.py','compare_sompy_features.py'):
                with self.subTest(script=script):
                    pd.DataFrame([row, dict(row, Filter='PASS')]).to_csv(a,index=False)
                    b.write_bytes(a.read_bytes())
                    command = [sys.executable,'-I','-B',str(ROOT/'src/sh'/script),str(a),str(b)]
                    result = subprocess.run(command,capture_output=True,text=True)
                    self.assertEqual(result.returncode,0,result.stderr)
                    changed=pd.read_csv(a)
                    changed.loc[1,'QUERY.TP']=4
                    changed.to_csv(b,index=False)
                    self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)
            metrics = ['fp','ambiguous','recall','recall_lower','recall_upper','precision','precision_lower','precision_upper','tp','total.query','ambi','na','recall2','unk','total.truth','fn','fp.region.size','fp.rate']
            record={'type':'SNVs', **dict.fromkeys(metrics,1)}
            pd.DataFrame([record]).to_csv(a)
            b.write_bytes(a.read_bytes())
            command=[sys.executable,'-I','-B',str(ROOT/'src/sh/compare_sompy.py'),str(a),str(b)]
            self.assertEqual(subprocess.run(command,capture_output=True).returncode,0)
            record['tp']=2
            pd.DataFrame([record]).to_csv(b)
            self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)

    def test_installed_layout_and_editor_helper(self):
        import runpy
        # ycm_core is an editor plugin, not a runtime dependency; stub only its
        # unused module while exercising the pure Python flag helper.
        with patch.dict(sys.modules, {'ycm_core': SimpleNamespace()}):
            config = runpy.run_path(str(ROOT/'.ycm_extra_conf.py'))
            flags = config['FlagsForFile']('example.cpp')['flags']
            self.assertTrue(flags)
            self.assertFalse(any('python2' in value for value in flags))
        with tempfile.TemporaryDirectory(prefix='installed 峰 ') as directory:
            prefix=Path(directory)
            bindir=prefix/'bin'
            modules=prefix/'lib/hap.py'
            bindir.mkdir()
            modules.mkdir(parents=True)
            for name in ('Tools','Haplo','Somatic'):
                shutil.copytree(ROOT/'src/python'/name,modules/name,ignore=shutil.ignore_patterns('__pycache__'))
            template=(modules/'Haplo/version.py.in').read_text()
            (modules/'Haplo/version.py').write_text(template.replace('@HAPLOTYPES_VERSION@','migration-test').replace('@SGE_AVAILABLE@','False').replace('@VCFEVAL_AVAILABLE@','False'))
            for name in ('hap.py','pre.py','qfy.py','som.py','ftx.py','cnx.py','bamstats.py','ovc.py'):
                shutil.copy2(ROOT/'src/python'/name,bindir/name)
            for name in ('hap.py','pre.py','qfy.py','som.py','ftx.py','cnx.py','bamstats.py','ovc.py'):
                result=subprocess.run([sys.executable,'-I','-B',str(bindir/name),'--help'],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
            for name in ('hap.py','pre.py'):
                result=subprocess.run([sys.executable,'-I','-B',str(bindir/name),'--version'],capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertIn('migration-test',result.stdout)

    def test_remaining_validation_scripts(self):
        result=subprocess.run([sys.executable,'-I','-B',str(ROOT/'src/sh/run_fastasize_test.py')],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        with tempfile.TemporaryDirectory() as directory:
            stats, happy = Path(directory)/'stats.csv', Path(directory)/'happy.csv'
            for script, stat_type, subset in [('validate_happy_summary.py','SNVs','*'),('validate_happy_extended.py','SNVs.0.00-0.50','[0.00,0.50)')]:
                row={'type':stat_type,'total.truth':3,'tp':2,'fn':1,'total.query':4,'fp':1,'unk':1,'ambi':0}
                pd.DataFrame([row]).to_csv(stats,index=False)
                output={'Type':'SNP','Subset':subset,'Filter':'ALL','TRUTH.TOTAL':3,'TRUTH.TP':2,'TRUTH.FN':1,'QUERY.TOTAL':4,'QUERY.FP':1,'QUERY.UNK':1}
                pd.DataFrame([output]).to_csv(happy,index=False)
                option='--happy_summary' if 'summary' in script else '--happy_extended'
                command=[sys.executable,'-I','-B',str(ROOT/'src/sh'/script),'--sompy_stats',str(stats),option,str(happy)]
                result=subprocess.run(command,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
                output['TRUTH.TP']=7
                pd.DataFrame([output]).to_csv(happy,index=False)
                self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)
                pd.DataFrame(columns=output.keys()).to_csv(happy,index=False)
                self.assertNotEqual(subprocess.run(command,capture_output=True).returncode,0)

    def test_vcf_independent_sample_fields(self):
        from Tools.vcfextract import vcfExtract
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'case.vcf.gz'
            with gzip.open(path, 'wt') as stream:
                stream.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\ncustom\t1\t.\tA\tC\t10\tPASS\t.\tGT:BD:BVT\t0/1:FN:SNP\t0/1:FP:SNP\n')
            self.assertEqual(list(vcfExtract(str(path), ['S.1.BD','S.2.BD','S.1.BVT','S.2.BVT'])), [['FN','FP','SNP','SNP']])


if __name__ == '__main__':
    from importlib.metadata import version
    print('Interpreter:', sys.executable, flush=True)
    print('Python:', sys.version, flush=True)
    for package in ('numpy', 'pandas', 'scipy', 'pysam', 'bx-python', 'psutil'):
        print(package, version(package), flush=True)
    print('pysam embedded HTSlib:', pysam.__samtools_version__, flush=True)
    unittest.main()
