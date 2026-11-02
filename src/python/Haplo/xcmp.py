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
# 3/9/2015
#
# Diploid VCF File Comparison
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import os
import logging
import tempfile
import time
import subprocess
import shlex


def xcmpWrapper(location_str, args):
    """ Haplotype block comparison wrapper function
    """
    starttime = time.time()
    tf = tempfile.NamedTemporaryFile(delete=False,
                                     dir=args.scratch_prefix,
                                     prefix="result.%s" % location_str,
                                     suffix=".bcf")
    tf.close()

    to_run = ['xcmp', args.vcf1, args.vcf2, '-l', location_str,
              '-o', tf.name, '-r', args.ref, '-f', str(int(args.pass_only)),
              '-n', str(args.max_enum), '--expand-hapblocks', str(args.hb_expand),
              '--window', str(args.window), '--no-hapcmp', str(int(args.no_hc)),
              '--qq', args.roc or 'QUAL']
    if args.verbose:
        to_run += ['-e', '-']

    # regions / targets already have been taken care of in blocksplit / preprocessing

    tfe = tempfile.NamedTemporaryFile(delete=False,
                                      dir=args.scratch_prefix,
                                      prefix="stderr",
                                      suffix=".log")
    tfo = tempfile.NamedTemporaryFile(delete=False,
                                      dir=args.scratch_prefix,
                                      prefix="stdout",
                                      suffix=".log")

    try:
        logging.info("Running '%s'" % to_run)
        subprocess.check_call(to_run, stdout=tfo, stderr=tfe)
    finally:
        tfo.close()
        tfe.close()
        with open(tfo.name) as f:
            for l in f:
                logging.info(l.replace("\n", ""))
        os.unlink(tfo.name)
        with open(tfe.name) as f:
            for l in f:
                logging.warning(l.replace("\n", ""))
        os.unlink(tfe.name)

    elapsed = time.time() - starttime
    logging.info("xcmp for chunk %s -- time taken %.2f" % (location_str, elapsed))

    return tf.name
