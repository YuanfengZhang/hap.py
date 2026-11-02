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
# 01/04/2015
#
# Process raw counts coming out of quantify

import os
import tempfile
import subprocess
import copy
import json
import logging
import Tools
import shlex

from Tools.bcftools import runBcftools


def _locations_tmp_bed_file(locations):
    """ turn a list of locations into a bed file """
    if type(locations) is str:
        locations = locations.split(",")
    if type(locations) is not list:
        raise Exception("Invalid list of locations (must be str or list): %s" % str(locations))

    llocations = []

    for l in locations:
        xchr, _, _pos = l.partition(":")
        start, _, end = _pos.partition("-")
        if not xchr:
            raise Exception("Invalid chromosome name in %s" % str(l))
        try:
            start = int(start)
        except:
            start = 0

        try:
            end = int(end)
        except:
            end = 2 ** 31 - 1

        llocations.append((xchr, start, end))

    locations = sorted(llocations)

    tf = tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", delete=False)
    for xchr, start, end in locations:
        print("%s\t%i\t%i" % (xchr, start - 1, end), file=tf)
    tf.close()

    return tf.name


def run_quantify(filename: str,
                 output_file=None, write_vcf=False, regions=None,
                 reference=None,
                 locations=None, threads=1,
                 output_vtc=False,
                 output_rocs=False,
                 qtype=None,
                 roc_file=None,
                 roc_val=None,
                 roc_header=None,
                 roc_filter=None,
                 roc_delta=None,
                 roc_regions=None,
                 clean_info=True,
                 strat_fixchr: bool = False) -> None:
    """Run native quantification and write the requested reports

    :param filename: the VCF file name
    :param output_file: output file name (if None, will use a temp file)
    :param write_vcf: write annotated VCF (give filename)
    :type write_vcf: str
    :param regions: dictionary of stratification region names and file names
    :param reference: reference fasta path
    :param locations: a location to use
    :param output_vtc: enable / disable the VTC field
    :param output_rocs: enable / disable output of ROCs by QQ level
    :param roc_file: filename for a TSV file with ROC observations
    :param roc_val: field to use for ROC QQ
    :param roc_header: name of ROC value for tables
    :param roc_filter: ROC filtering settings
    :param roc_delta: ROC minimum spacing between levels
    :param roc_regions: List of regions to output full ROCs for
    :param clean_info: remove unused INFO fields
    :param strat_fixchr: fix chr naming in stratification regions
    :returns: None
    """

    reference = Tools.defaultReference(reference)
    if reference is None:
        raise ValueError("Quantification requires a reference FASTA (-r or HGREF).")

    if not output_file:
        output_file = tempfile.NamedTemporaryFile().name

    command: list[str] = ['quantify', filename, '-o', output_file,
                          '-r', reference, '--threads', str(threads),
                          '--output-vtc', str(int(bool(output_vtc))),
                          '--output-rocs', str(int(bool(output_rocs))),
                          '--clean-info', str(int(bool(clean_info))),
                          '--fix-chr-regions', str(int(bool(strat_fixchr)))]
    for option, value in (('--type', qtype), ('--output-roc', roc_file),
                          ('--qq', roc_val), ('--roc-filter', roc_filter)):
        if value:
            command.extend((option, value))
    if roc_val and roc_header != roc_val:
        command.extend(('--qq-header', roc_header))
    if roc_delta:
        command.extend(('--roc-delta', '%f' % roc_delta))
    if write_vcf:
        if not write_vcf.endswith(('.vcf.gz', '.bcf')):
            write_vcf += '.vcf.gz'
        command.extend(('-v', write_vcf))
    for name, path in (regions or {}).items():
        command.extend(('-R', name + ':' + path))
    for region in roc_regions or ():
        command.extend(('--roc-regions', region))
    location_file = None
    if locations:
        location_file = _locations_tmp_bed_file(locations)
        command.extend(('--only', location_file))

    tfe = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stderr",
                                      suffix=".log")
    tfo = tempfile.NamedTemporaryFile(delete=False,
                                      prefix="stdout",
                                      suffix=".log")

    logging.info('Running %s', shlex.join(command))

    try:
        subprocess.check_call(command, stdout=tfo, stderr=tfe)
    except:
        tfo.close()
        tfe.close()
        with open(tfo.name) as f:
            for l in f:
                logging.error("[stdout] " + l.replace("\n", ""))
        os.unlink(tfo.name)
        with open(tfe.name) as f:
            for l in f:
                logging.error("[stderr] " + l.replace("\n", ""))
        os.unlink(tfe.name)
        if location_file:
            os.unlink(location_file)
        raise

    tfo.close()
    tfe.close()
    with open(tfo.name) as f:
        for l in f:
            logging.info("[stdout] " + l.replace("\n", ""))
    os.unlink(tfo.name)
    with open(tfe.name) as f:
        for l in f:
            logging.info("[stderr] " + l.replace("\n", ""))
    os.unlink(tfe.name)
    if location_file:
        os.unlink(location_file)

    if write_vcf and write_vcf.endswith(".bcf"):
        runBcftools("index", write_vcf)
    elif write_vcf:
        command = ['tabix', '-p', 'vcf', write_vcf]
        logging.info('Running %s', shlex.join(command))
        subprocess.check_call(command)


