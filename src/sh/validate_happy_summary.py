#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Compare som.py stats.csv and hap.py summary.csv files

import sys
import argparse
import csv
import pprint as pp
import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO)

def parse_args():
    parser = argparse.ArgumentParser("Validate hap.py stats")
    parser.add_argument("--sompy_stats", required=True, help="Path to som.py stats.csv")
    parser.add_argument("--happy_summary", required=True, help="Path to hap.py summary.csv")
    args = parser.parse_args()
    return args

def eval_equal(metric_name, count_a, count_b):
    a = int(count_a)
    b = int(count_b)
    e = "PASS" if a == b else "FAIL"
    logging.info("%s: %d vs %d - %s" % (metric_name, a, b, e))
    return e

if __name__ == '__main__':
    args = parse_args()
    sompy_stats = csv.DictReader(open(args.sompy_stats))
    happy_summary = csv.DictReader(open(args.happy_summary))

    # records counts VCF rows, not a disjoint variant class (bcftools SN).
    # Validate every typed row; require a common matching ALL or PASS filter.
    happy_summary = list(happy_summary)
    outcomes = dict(ALL=set(), PASS=set())
    seen_types = set()
    for s in sompy_stats:
        if s["type"] == "records":
            continue
        if s["type"] == "SNVs":
            vtype = "SNP"
        elif s["type"] == "indels":
            vtype = "INDEL"
        else:
            logging.error("Unrecognised variant type in som.py stats: %s" % s["type"])
            sys.exit(1)

        seen_types.add(vtype)
        logging.info("Comparing %s counts..." % vtype)
        matched_filters = set()
        for h in happy_summary:
            if h["Type"] == vtype and h["Filter"] in outcomes:
                matched_filters.add(h["Filter"])
                outcomes[h["Filter"]].add(eval_equal(metric_name="TRUTH.TOTAL", count_a=s["total.truth"], count_b=h["TRUTH.TOTAL"]))
                outcomes[h["Filter"]].add(eval_equal(metric_name="TRUTH.TP", count_a=s["tp"], count_b=h["TRUTH.TP"]))
                outcomes[h["Filter"]].add(eval_equal(metric_name="TRUTH.FN", count_a=s["fn"], count_b=h["TRUTH.FN"]))
                outcomes[h["Filter"]].add(eval_equal(metric_name="QUERY.TOTAL", count_a=s["total.query"], count_b=h["QUERY.TOTAL"]))
                outcomes[h["Filter"]].add(eval_equal(metric_name="QUERY.FP", count_a=s["fp"], count_b=h["QUERY.FP"]))
                outcomes[h["Filter"]].add(eval_equal(metric_name="QUERY.UNK", count_a=int(s["unk"])+int(s["ambi"]), count_b=h["QUERY.UNK"]))

        for vfilter in outcomes:
            if vfilter not in matched_filters:
                outcomes[vfilter].add("FAIL")

    if not seen_types:
        logging.error("No SNVs or indels rows in som.py stats")
        sys.exit(1)

    failed_vfilters = [x for x in outcomes if "FAIL" in outcomes[x]]
    if not any(values == {"PASS"} for values in outcomes.values()):
        logging.info("Failed filters: %s" % failed_vfilters)
        sys.exit(1)
    else:
        logging.info("DONE")
