#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This is a modified version of the original software.
# The original authors' notice follows.
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# 10/1/2015
#
# Testing helper, will return exit code 1 if the passed bed file has
# overlapping intervals.
#
# Usage:
#
# ovc.py input.bed
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import argparse
import gzip


def main() -> int:
    parser = argparse.ArgumentParser(description='Check BED intervals for overlaps')
    parser.add_argument('input')
    args = parser.parse_args()
    opener = gzip.open if args.input.endswith('.gz') else open
    last, chromosome = -1, None
    with opener(args.input, 'rt', encoding='utf-8') as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                raise ValueError(f'Invalid BED at line {number}')
            start, end = int(fields[1]), int(fields[2])
            if fields[0] != chromosome:
                last, chromosome = -1, fields[0]
            if start < last:
                print(f'Overlap at {chromosome}:{start} (line {number})')
                return 1
            last = end
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
