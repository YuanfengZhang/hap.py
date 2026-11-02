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
import errno
import logging
import subprocess
from datetime import date

# noinspection PyUnresolvedReferences
try:
    import Haplo.version as vs  # pylint: disable=E0401,E0611
    version = vs.__version__
    has_sge = vs.has_sge
    has_vcfeval = vs.has_vcfeval
except ImportError:
    logging.warning("No version found. Please follow the installation instructions.")
    version = "unknown"
    has_sge = False
    has_vcfeval = False


def defaultReference(reference: str | None = None) -> str | None:
    """Resolve explicit reference > HGREF > HG19 > conventional data path.

    A configured but invalid reference is an error, never a silent fallback to
    a different assembly. With no configuration, absence returns None so --help
    remains usable without reference data.
    """
    source = 'explicit reference'
    candidate = reference
    if not candidate:
        for key in ('HGREF', 'HG19'):
            if os.environ.get(key):
                source, candidate = key, os.environ[key]
                break
    if not candidate:
        candidate = '/opt/hap.py-data/hg19.fa'
        if not os.path.isfile(candidate):
            return None
    candidate = os.path.abspath(os.path.expanduser(candidate))
    if not os.path.isfile(candidate) or not os.access(candidate, os.R_OK):
        raise ValueError(f'{source} FASTA is not readable: {candidate}. '
                         'Use -r/--reference or set HGREF to the correct FASTA.')
    return candidate


def which(program):
    def is_exe(xfpath):
        return os.path.isfile(xfpath) and os.access(xfpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def init():
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)

    base = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
    paths = ["bin"]

    for p in paths:
        pp = os.path.join(base, p)
        if not os.path.exists(pp):
            raise Exception("Dependency path %s not found" % pp)
        os.environ["PATH"] = pp + os.pathsep + os.environ["PATH"]

    # Validate tools when invoked: each CLI uses a different subset, and
    # importing helpers or running --help must not require unused executables.

    os.environ['DYLD_LIBRARY_PATH'] = os.path.join(base, "lib")
    os.environ['LD_LIBRARY_PATH'] = os.path.join(base, "lib")


# Initialization is explicit at CLI execution, never during worker imports.

# safely import here
# noinspection PyUnresolvedReferences
import pysam  # noqa: F401
# noinspection PyUnresolvedReferences
import pandas  # noqa: F401


class LoggingWriter(object):
    """ Helper class to write tracebacks to log file
    """
    def __init__(self, level):
        self.level = level

    def write(self, message):
        message = message.replace("\n", "")
        if message:
            logging.log(self.level, message)


def writeVCFHeader(filelike, extrainfo='', chrprefix=None, reference: str | None = None) -> None:
    """Write a VCF header using the actual FASTA sequence dictionary.

    chrprefix is retained for call compatibility but may not rename reference
    sequences. FASTA contig identifiers are authoritative.
    """
    reference = defaultReference(reference)
    if reference is None:
        raise ValueError('VCF header requires a FASTA: pass reference or set HGREF.')
    with pysam.FastaFile(reference) as fasta:
        contigs = list(zip(fasta.references, fasta.lengths))
    extra = extrainfo if isinstance(extrainfo, list) else extrainfo.splitlines()
    if any(line.startswith(('##reference=', '##contig=')) for line in extra):
        raise ValueError('Reference and contig metadata must come from the FASTA, not extrainfo.')
    infos = ['##fileformat=VCFv4.2', f'##reference={reference}',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    infos.extend(extra)
    infos.extend(f'##contig=<ID={name},length={length}>' for name, length in contigs)
    infos.extend(['##fileDate=' + date.today().isoformat(),
                  '##source=HaploCompare', f'##source_version={version}'])
    for line in infos:
        filelike.write(line + '\n')
    filelike.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIMPLE\n')


def mkdir_p(path):
    """ mkdir -p path """
    try:
        os.makedirs(os.path.abspath(path))
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
    if not os.path.isdir(path):
        raise Exception("Failed to create directory %s" % path)


class BGZipFile:
    """BGZF writer accepting UTF-8 text or raw bytes, without shell pipes."""

    def __init__(self, filename: str, force: bool = False) -> None:
        if os.path.exists(filename) and not force:
            raise FileExistsError(f'File {filename} exists, use force=True to overwrite')
        self.name = filename
        self.write_file = pysam.BGZFile(filename, 'wb')

    def close(self) -> None:
        self.write_file.close()

    def write(self, data: str | bytes) -> int:
        return self.write_file.write(data.encode('utf-8') if isinstance(data, str) else data)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()
