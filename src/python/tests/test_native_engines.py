#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Run actual installed hap.py engines against a minimal reference fixture.

Usage: python -I -B test_native_engines.py NATIVE_BIN [RTG_EXECUTABLE]
Requires bcftools/samtools/tabix on PATH. No fake native executables are used.
"""
import json
import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'src/python'))
import pandas as pd
import pysam


def main() -> None:
    native = Path(sys.argv[1]).resolve()
    start_method = os.environ.get('HAPPY_TEST_START_METHOD')
    roc_enabled = os.environ.get('HAPPY_TEST_ROC') == '1'
    stage = ROOT/'.work/installed-integration'
    stage.mkdir(exist_ok=True)
    bindir, modules = stage/'bin', stage/'lib/hap.py'
    bindir.mkdir(exist_ok=True)
    modules.mkdir(parents=True, exist_ok=True)
    for name in ('Tools', 'Haplo', 'Somatic'):
        shutil.copytree(ROOT/'src/python'/name, modules/name, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__'))
    template = (modules/'Haplo/version.py.in').read_text()
    (modules/'Haplo/version.py').write_text(template.replace('@HAPLOTYPES_VERSION@','migration-test').replace('@SGE_AVAILABLE@','False').replace('@VCFEVAL_AVAILABLE@','False'))
    for path in (ROOT/'src/python').glob('*.py'):
        shutil.copy2(path, bindir/path.name)
    environment = dict(os.environ)
    environment['PATH'] = str(native)+os.pathsep+environment['PATH']
    fixture = stage/('fixture roc unk 峰' if roc_enabled else 'fixture space 峰')
    if start_method:
        fixture = fixture.with_name(fixture.name+' '+start_method)
    fixture.mkdir(exist_ok=True)
    reference=fixture/'reference.fa'
    reference.write_text('>custom\n'+'ACGT'*100+'\n')
    pysam.faidx(str(reference))
    header='##fileformat=VCFv4.2\n##contig=<ID=custom,length=400>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    truth=fixture/'truth.vcf'
    query=fixture/'query.vcf'
    truth.write_text(header+'custom\t2\t.\tC\tT\t40\tPASS\t.\tGT\t0/1\ncustom\t10\t.\tC\tA\t20\tPASS\t.\tGT\t0/1\n')
    query.write_text(header+'custom\t2\t.\tC\tT\t40\tPASS\t.\tGT\t0/1\ncustom\t18\t.\tC\tA\t10\tPASS\t.\tGT\t0/1\n')
    if os.environ.get('HAPPY_TEST_PRESERVE_INFO') == '1':
        for source, value in ((truth, 7), (query, 11)):
            content = source.read_text().replace('##FORMAT=', '##INFO=<ID=SOURCE_SCORE,Number=1,Type=Integer,Description="Source score">\n##FORMAT=', 1)
            source.write_text(content.replace('\tPASS\t.\tGT', f'\tPASS\tSOURCE_SCORE={value}\tGT'))
    if roc_enabled:
        with query.open('a') as stream:
            stream.write('custom\t302\t.\tC\tA\t2\tPASS\t.\tGT\t0/1\n')
    bed=fixture/'confident.bed'
    bed.write_text('custom\t0\t100\n' if roc_enabled else 'custom\t0\t400\n')
    engines=os.environ.get('HAPPY_TEST_ENGINES','xcmp,scmp-somatic').split(',')
    if len(sys.argv)>2:
        engines.append('vcfeval')
    for engine in engines:
        for threads in (1,2):
            prefix=fixture/f'{engine}-{threads}'
            if roc_enabled:
                prefix = fixture/f'{engine}-{threads}-thresholds'
            command=[sys.executable,'-I','-B',str(bindir/'hap.py'),str(truth),str(query),'-r',str(reference),'-f',str(bed),'-o',str(prefix),'--engine',engine,'--threads',str(threads),'--scratch-prefix',str(fixture),'--no-roc','--verbose','--keep-scratch']
            if roc_enabled:
                command.remove('--no-roc')
                command.append('-V')
            if os.environ.get('HAPPY_TEST_PRESERVE_INFO') == '1':
                command.append('--preserve-info')
            if engine=='vcfeval':
                command += ['--engine-vcfeval-path',sys.argv[2]]
            if start_method:
                command = [sys.executable, '-I', '-B', '-c',
                           'import multiprocessing,runpy,sys; multiprocessing.set_start_method(sys.argv.pop(1)); sys.argv.pop(0); runpy.run_path(sys.argv[0], run_name="__main__")',
                           start_method, *command[3:]]
            log=ROOT/'.work/logs'/f'engine-{engine}-{threads}{"-thresholds" if roc_enabled else ""}{"-"+start_method if start_method else ""}.log'
            print(json.dumps(command,ensure_ascii=False),flush=True)
            with log.open('w') as stream:
                result=subprocess.run(command,env=environment,stdout=stream,stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError(f'{engine}/{threads} exited {result.returncode}; see {log}')
            summary=pd.read_csv(str(prefix)+'.summary.csv')
            row=summary[(summary['Type']=='SNP') & (summary['Filter']=='ALL')].iloc[0]
            for field,expected in {'TRUTH.TOTAL':2,'TRUTH.TP':1,'TRUTH.FN':1,'QUERY.TOTAL':3 if roc_enabled else 2,'QUERY.FP':1,'QUERY.UNK':1 if roc_enabled else 0,'METRIC.Recall':.5,'METRIC.Precision':.5,'METRIC.F1_Score':.5}.items():
                if row[field] != expected:
                    raise AssertionError(f'{engine}/{threads} {field}: {row[field]} != {expected}')
            print(engine,threads,'exact counts and metrics PASS',flush=True)
            if os.environ.get('HAPPY_TEST_PRESERVE_INFO') == '1':
                with pysam.VariantFile(str(prefix)+'.vcf.gz') as records:
                    if tuple(records.header.samples) != ('TRUTH', 'QUERY'):
                        raise AssertionError('INFO preservation changed comparison samples')
                    annotated = [record for record in records if record.pos == 2]
                    if len(annotated) != 1 or annotated[0].info.get('SOURCE_SCORE') != 7:
                        raise AssertionError('Original INFO was not preserved on the shared variant')

            if roc_enabled:
                table = pd.read_csv(str(prefix)+'.roc.all.csv.gz')
                for _, record in table.iterrows():
                    truth_den = record['TRUTH.TP'] + record['TRUTH.FN']
                    query_den = record['QUERY.TP'] + record['QUERY.FP']
                    recall = record['TRUTH.TP']/truth_den if truth_den else 0.0
                    precision = record['QUERY.TP']/query_den if query_den else 0.0
                    if not (abs(record['METRIC.Recall']-recall) < 1e-12):
                        raise AssertionError("Regression condition failed: abs(record['METRIC.Recall']-recall) < 1e-12")
                    if not (abs(record['METRIC.Precision']-precision) < 1e-12):
                        raise AssertionError("Regression condition failed: abs(record['METRIC.Precision']-precision) < 1e-12")
                empty = summary[summary['Type']=='INDEL']
                # Absent variant types are omitted by the native output contract.
                if empty.empty:
                    print(engine, 'empty INDEL stratum omitted (no synthetic row)', flush=True)
                if not ((empty[['TRUTH.TP','TRUTH.FN','QUERY.FP','METRIC.Recall','METRIC.Precision','METRIC.F1_Score']]==0).all().all()):
                    raise AssertionError("Regression condition failed: (empty[['TRUTH.TP','TRUTH.FN','QUERY.FP','METRIC.Recall','METRIC.Precision','METRIC.F1_Score']]==0).all().all()")
                selected = table[(table['Type']=='SNP') & (table['Subtype']=='*') & (table['Genotype']=='*') & (table['Subset']=='*') & (table['Filter']=='ALL')]
                thresholds = selected[selected['QQ']!='*'].copy()
                numeric = pd.to_numeric(thresholds['QQ'])
                print(engine, threads, 'thresholds', numeric.tolist(), flush=True)
                if not (set(numeric)=={0.0,2.0,10.0,40.0}):
                    raise AssertionError('Regression condition failed: set(numeric)=={0.0,2.0,10.0,40.0}')
                if not (numeric.is_monotonic_increasing):
                    raise AssertionError('ROC thresholds not numerically sorted')
                thresholds.index = numeric
                for threshold, expected_fp, expected_unk in ((0,1,1),(2,1,1),(10,1,0),(40,0,0)):
                    record = thresholds.loc[threshold]
                    for field, expected in {'TRUTH.TP':1,'TRUTH.FN':1,'QUERY.TP':1,'QUERY.FP':expected_fp,'QUERY.UNK':expected_unk}.items():
                        if not (record[field] == expected):
                            raise AssertionError((engine, threshold, field, record[field], expected))
                import gzip
                with gzip.open(str(prefix)+'.metrics.json.gz', 'rt') as stream:
                    metrics = json.load(stream)
                for metric in metrics['metrics']:
                    for column in metric.get('data', []):
                        if column['id'] in ('Subset.Size','Subset.IS_CONF.Size','Subset.Level'):
                            if not (column['type']=='int64'):
                                raise AssertionError("Regression condition failed: column['type']=='int64'")
                            if not (all(value is None or type(value) is int for value in column['values'])):
                                raise AssertionError("Regression condition failed: all(value is None or type(value) is int for value in column['values'])")
                print(engine,threads,'thresholds, empty strata and independent denominators PASS',flush=True)
        suffix = '-thresholds' if roc_enabled else ''
        left=pd.read_csv(str(fixture/f'{engine}-1{suffix}')+'.extended.csv')
        right=pd.read_csv(str(fixture/f'{engine}-2{suffix}')+'.extended.csv')
        pd.testing.assert_frame_equal(left,right)
        if roc_enabled:
            import gzip
            prefixes = [fixture/f'{engine}-{number}{suffix}' for number in (1,2)]
            def read_vcf(prefix: Path) -> tuple:
                with pysam.VariantFile(str(prefix)+'.vcf.gz') as stream:
                    contigs = tuple((name, item.length) for name, item in stream.header.contigs.items())
                    definitions = tuple(str(item) for item in stream.header.records if item.key in ('INFO', 'FORMAT', 'FILTER'))
                    records = tuple((record.contig, record.pos, record.stop, record.id, record.ref, record.alts,
                                  record.qual, tuple(record.filter), dict(record.info),
                                  {name: (dict(sample), sample.phased) for name, sample in record.samples.items()})
                                 for record in stream)
                    return contigs, definitions, tuple(stream.header.samples), records
            if not (read_vcf(prefixes[0]) == read_vcf(prefixes[1])):
                raise AssertionError('Regression condition failed: read_vcf(prefixes[0]) == read_vcf(prefixes[1])')
            payloads = []
            for prefix in prefixes:
                with gzip.open(str(prefix)+'.metrics.json.gz', 'rt') as stream:
                    payload = json.load(stream)
                payload.pop('timestamp', None)
                payload['runInfo'] = [item for item in payload.get('runInfo', []) if item.get('key') != 'commandline']
                required = payload.get('metadata', {}).get('required', {})
                if 'description' in required:
                    required['description'] = required['description'].split(' generated this JSON file via command line ')[0]
                payloads.append(payload)
            if not (payloads[0] == payloads[1]):
                raise AssertionError('Regression condition failed: payloads[0] == payloads[1]')
            outputs = list(fixture.glob(prefixes[0].name+'.roc.*.csv.gz'))
            if not (outputs):
                raise AssertionError('Regression condition failed: outputs')
            for first in outputs:
                second = fixture/first.name.replace(prefixes[0].name, prefixes[1].name, 1)
                pd.testing.assert_frame_equal(pd.read_csv(first), pd.read_csv(second))
            print(engine, 'VCF records, JSON metrics and all ROC tables thread-equivalent PASS', flush=True)
        print(engine,'single/multiprocess equivalence PASS',flush=True)
        if start_method:
            failed = list(command)
            failed[failed.index('-r')+1] = str(fixture/'missing-reference.fa')
            failure_log = ROOT/'.work/logs'/f'engine-{engine}-{start_method}-failure.log'
            with failure_log.open('w') as stream:
                result = subprocess.run(failed, env=environment, stdout=stream, stderr=subprocess.STDOUT)
            if not (result.returncode != 0):
                raise AssertionError('CLI must propagate invalid reference failure')
            print(engine, start_method, 'CLI invalid-reference failure PASS', flush=True)


if __name__=='__main__':
    main()
