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

import os
import sys
import gc
import argparse
import logging
import traceback
import subprocess
import multiprocessing
import pickle
import tempfile
from itertools import islice, repeat

from . import LoggingWriter


def getPool(threads: int):
    """Return an independently owned pool; callers must close/join it."""
    return multiprocessing.get_context().Pool(threads) if threads > 1 else None


def splitEvery(n: int | None, iterable):
    """Yield finite blocks, rejecting nonpositive block sizes."""
    if n is None:
        yield list(iterable)
        return
    if n <= 0:
        raise ValueError('Block size must be positive')
    iterator = iter(iterable)
    while piece := list(islice(iterator, n)):
        yield piece


def unpickleSequentially(plist):
    """ Unpickle and concatenate sequentially """
    data = []
    while plist or data:
        if not data:
            fname = plist.pop(0)
            with open(fname, 'rb') as f:
                data = pickle.load(f)
            os.unlink(fname)
        if data:
            yield data.pop(0)

def parMapper(arg):
    """Let worker failures propagate to the parent instead of returning None."""
    gc.collect()
    return arg[1]['fun'](arg[0], *arg[1]['args'], **arg[1]['kwargs'])


def runParallel(pool, fun, par, *args, **kwargs):
    """ run a function in parallel on all elements in par

    :param pool: multiprocessing.Pool or None
    :param fun: a function
    :param par: a list of things to map to (each item is passed as the first argument to fun)
    :param args: more function arguments for fun
    :param kwargs: more function arguments for fun

    """
    if pool:
        result = pool.map(parMapper, zip(par, repeat( { "fun": fun, "args": args, "kwargs": kwargs } )))
    else:
        result = []
        for c in par:
            result.append(parMapper( (c, { "fun": fun, "args": args, "kwargs": kwargs } ) ))
    return result
