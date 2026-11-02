#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare all summary strata, including independent truth/query counts."""
import argparse
import csv
from decimal import Decimal, InvalidOperation
from pathlib import Path


def csvread(filename: str) -> dict:
    data = {}
    with open(filename, newline='', encoding='utf-8') as stream:
        for record in csv.DictReader(stream):
            first = next(iter(record))
            if 'hap.py' in record[first]:
                continue
            label = tuple(record.get(key, '') for key in
                          (first, 'Subtype', 'Subset', 'Filter', 'Genotype', 'QQ.Field', 'QQ'))
            if label in data:
                raise ValueError(f'Duplicate summary stratum: {label}')
            data[label] = record
    return data


def validate_schema(record: dict, expected: dict, extended: list[dict], qq_field: str) -> None:
    added: set[str] = record.keys() - expected.keys()
    if expected.keys() - record.keys() or not added <= {'FP.al', 'METRIC.F1_Score'}:
        raise AssertionError('Unrecognized summary schema change')
    matches: list[dict] = [row for row in extended
                          if row['Type'] == record['Type'] and row['Filter'] == record['Filter']
                          and row.get('QQ.Field') == qq_field
                          and all(row.get(key) == '*' for key in ('Subtype', 'Subset', 'Genotype', 'QQ'))]
    if len(matches) != 1:
        raise AssertionError('Expected one matching extended stratum')
    row: dict = matches[0]
    for column, value in record.items():
        if column.startswith(('TRUTH.', 'QUERY.', 'FP.')) and not column.endswith('_ratio'):
            if Decimal(value) != Decimal(row[column]):
                raise AssertionError(f'Summary/extended count differs: {column}')
    tp, fn, query_tp, fp = (Decimal(row[key]) for key in
                            ('TRUTH.TP', 'TRUTH.FN', 'QUERY.TP', 'QUERY.FP'))
    recall: Decimal = tp / (tp + fn) if tp + fn else Decimal(0)
    precision: Decimal = query_tp / (query_tp + fp) if query_tp + fp else Decimal(0)
    f1: Decimal = 2 * precision * recall / (precision + recall) if precision + recall else Decimal(0)
    if abs(Decimal(record['METRIC.F1_Score']) - f1) > Decimal('1e-14'):
        raise AssertionError('Summary F1 differs from independent truth/query denominators')
    allele_fp: Decimal = Decimal(record['FP.al'])
    if not allele_fp.is_finite() or allele_fp != allele_fp.to_integral_value() or not 0 <= allele_fp <= fp:
        raise AssertionError('Invalid allele mismatch count')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('actual')
    parser.add_argument('expected')
    parser.add_argument('--schema', choices=['ga4gh-summary-v2'])
    parser.add_argument('--qq-field', default='QUAL', help='Expected extended-table score field')
    args = parser.parse_args()
    left, right = csvread(args.actual), csvread(args.expected)
    extended: list[dict] = []
    if args.schema:
        if not args.actual.endswith('.summary.csv'):
            parser.error('Versioned schema requires a .summary.csv and companion .extended.csv')
        companion: Path = Path(args.actual[:-len('.summary.csv')] + '.extended.csv')
        with companion.open(newline='', encoding='utf-8') as stream:
            extended = list(csv.DictReader(stream))
    if left.keys() != right.keys() or not left:
        raise AssertionError('Summary strata differ or are empty')
    for label, record in left.items():
        if args.schema:
            validate_schema(record, right[label], extended, args.qq_field)
        elif record.keys() != right[label].keys():
            raise AssertionError(f'Summary columns differ at {label}')
        for column in right[label]:
            value = record[column]
            other = right[label][column]
            if value == other:
                continue
            try:
                a, b = Decimal(value), Decimal(other)
            except InvalidOperation:
                if column.endswith('_ratio') and value in ('', '.') and other in ('', '.'):
                    continue
                raise AssertionError(f'{label} {column}: {value} != {other}')
            is_count = column.startswith(('TRUTH.', 'QUERY.', 'FP.')) and not column.endswith('_ratio')
            equal = a == b if is_count else ('%.3g' % a) == ('%.3g' % b)
            if not equal and not (a.is_nan() and b.is_nan()):
                raise AssertionError(f'{label} {column}: {value} != {other}')


if __name__ == '__main__':
    main()
