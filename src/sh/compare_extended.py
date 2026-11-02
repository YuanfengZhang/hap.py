#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Compare two extended csv files

import sys
import argparse
import csv
from decimal import Decimal

def csvread(filename):
    f = open(filename)
    freader = csv.DictReader(f)
    data = {}
    label_columns = ["Type", 
                     "Subtype", 
                     "Subset",
                     "Filter", 
                     "Genotype",
                     "QQ.Field",
                     "QQ"]

    for l in freader:
        record = dict(l)
        label = tuple(record[k] for k in label_columns)
        if "hap.py" in label:
            # ignore version line
            continue

        for k in list(record.keys()):
            try:
                record[k] = int(record[k])
            except:
                try:
                    record[k] = Decimal(record[k])
                except:
                    if record[k] == "." or record[k] == "":
                        record[k] = Decimal("NaN")
        if label in data:
            raise ValueError("Duplicate extended stratum: %s" % (label,))
        data[label] = record

    return data


def validate_metric_schema(actual: dict, expected: dict) -> None:
    tp, fn, qtp, fp, total, unk = (Decimal(actual[key]) for key in
        ('TRUTH.TP', 'TRUTH.FN', 'QUERY.TP', 'QUERY.FP', 'QUERY.TOTAL', 'QUERY.UNK'))
    if not all(value.is_finite() and value >= 0 for value in (tp, fn, qtp, fp, total, unk)):
        raise ValueError('Non-finite or negative metric count')
    recall = tp / (tp + fn) if tp + fn else Decimal(0)
    precision = qtp / (qtp + fp) if qtp + fp else Decimal(0)
    definitions = {
        'METRIC.Recall': (recall, tp + fn),
        'METRIC.Precision': (precision, qtp + fp),
        'METRIC.F1_Score': (2 * precision * recall / (precision + recall)
                            if precision + recall else Decimal(0), precision + recall),
        'METRIC.Frac_NA': (unk / total if total else Decimal(0), total),
    }
    for column, (value, denominator) in definitions.items():
        measured = Decimal(actual[column])
        if not measured.is_finite() or abs(measured - value) > Decimal('1e-14'):
            raise ValueError('Invalid metric formula: %s' % column)
        # Adapt only legacy missing values with independent zero-denominator evidence; compare all counts strictly below.
        if Decimal(expected[column]).is_nan() and denominator == 0:
            expected[column] = Decimal(0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('actual')
    parser.add_argument('expected')
    parser.add_argument('--schema', choices=['ga4gh-metrics-v2'])
    args = parser.parse_args()
    data1 = csvread(args.actual)
    data2 = csvread(args.expected)

    all_keys_1 = set()
    for l in list(data1.values()):
        all_keys_1 |= set(l.keys())
    all_keys_2 = set()
    for l in list(data2.values()):
        all_keys_2 |= set(l.keys())

    if len(all_keys_1 - all_keys_2) or len(all_keys_2 - all_keys_1):
        print("ERROR -- missing metrics:", file=sys.stderr)
        print("-------------------------\n", file=sys.stderr)
        print("In 1 but not 2: %s" % str(all_keys_1 - all_keys_2), file=sys.stderr)
        print("In 2 but not 1: %s" % str(all_keys_2 - all_keys_1), file=sys.stderr)
        print("-------------------------", file=sys.stderr)
        sys.exit(1)

    all_labels_1 = set(data1.keys())
    all_labels_2 = set(data2.keys())

    if len(all_labels_1 - all_labels_2) or len(all_labels_2 - all_labels_1):
        print("ERROR -- missing rows:", file=sys.stderr)
        print("-------------------------\n", file=sys.stderr)
        print("In 1 but not 2: %s" % str(all_labels_1 - all_labels_2), file=sys.stderr)
        print("In 2 but not 1: %s" % str(all_labels_2 - all_labels_1), file=sys.stderr)
        print("-------------------------", file=sys.stderr)
        sys.exit(1)

    print("Comparing %i labels and %i metrics..." % (len(all_labels_1), len(all_keys_1)))

    different_metrics = []

    for label in all_labels_1:
        if args.schema:
            validate_metric_schema(data1[label], data2[label])
        for metric in all_keys_1:
            v1 = data1[label][metric]
            v2 = data2[label][metric]

            if metric.endswith("_ratio"):
                if v1 == ".":
                    v1 = ""
                if v2 == ".":
                    v2 = ""

            if isinstance(v1, (int, Decimal)) and isinstance(v2, (int, Decimal)):
                a, b = Decimal(v1), Decimal(v2)
                exact = metric.startswith(('TRUTH.', 'QUERY.')) and not metric.endswith('_ratio')
                exact = exact or metric.startswith('FP.') or metric.startswith('Subset.')
                equal = (a.is_nan() and b.is_nan()) or a == b
                if not exact and not equal:
                    equal = ('%.3g' % a) == ('%.3g' % b)
                if not equal:
                    different_metrics.append((label, metric, str(v1), str(v2), str(a - b)))
            elif type(v1) != type(v2):
                different_metrics.append((label,
                                          metric,
                                          ("%s" % type(v1)),
                                          ("%s" % type(v2)),
                                          "type mismatch"))
            elif type(v1) is float:
                if ("%.3g" % v1) != ("%.3g" % v2):
                    different_metrics.append((label,
                                              metric,
                                              ("%.3g" % v1),
                                              ("%.3g" % v2),
                                              str(v1 - v2)))
            elif type(v1) is int:
                if v1 != v2:
                    different_metrics.append((label,
                                              metric,
                                              str(v1),
                                              str(v2),
                                              str(v1 - v2)))
            else:
                if v1 != v2:
                    different_metrics.append((label,
                                              metric,
                                              str(v1),
                                              str(v2),
                                              "Non-numeric mismatch"))

    if different_metrics:
        print("ERROR -- Metric differences detected:", file=sys.stderr)
        print("-------------------------------------\n", file=sys.stderr)
        for m in different_metrics:
            print("%s / %s: %s != %s difference: %s" % m, file=sys.stderr)
        print("-------------------------------------", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
