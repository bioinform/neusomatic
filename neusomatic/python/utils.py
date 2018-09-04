#-------------------------------------------------------------------------
# utils.py
# Utility functions
#-------------------------------------------------------------------------
import os
import shutil
import shlex
import subprocess
import logging

import pysam
import numpy as np


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)
logging.getLogger().setLevel(logging.INFO)


def run_shell_command(command, stdout=None, stderr=None):
    stdout_fd = open(stdout, "w") if stdout else None
    stderr_fd = open(stderr, "w") if stderr else None

    fixed_command = shlex.split(command)
    logger.info("Running command: {}".format(fixed_command))
    process = subprocess.Popen(
        fixed_command, stdout=stdout_fd, stderr=stderr_fd)
    process.wait()
    if stdout_fd:
        stdout_fd.close()
    if stderr_fd:
        stderr_fd.close()
    return process.returncode


def concatenate_files(infiles, outfile, check_file_existence=True):
    with open(outfile, "w") as out_fd:
        for infile in infiles:
            if not infile or check_file_existence and not os.path.isfile(infile):
                continue
            with open(infile) as in_fd:
                shutil.copyfileobj(in_fd, out_fd)
    return outfile


def concatenate_vcfs(infiles, outfile, check_file_existence=True, header_string="#"):
    with open(outfile, "w") as out_fd:
        # Only keep files which exist
        files_to_process = filter(lambda f: f and (
            not check_file_existence or os.path.isfile(f)), infiles)

        for index, infile in enumerate(files_to_process):
            with open(infile) as in_fd:
                if index == 0:
                    shutil.copyfileobj(in_fd, out_fd)
                else:
                    for line in in_fd:
                        if not line.startswith(header_string):
                            out_fd.write(line)
    return outfile


def get_chromosomes_order(reference=None, bam=None):
    chroms_order = {}
    if reference:
        with pysam.FastaFile(reference) as fd:
            chroms_order = {chrom: chrom_index for chrom_index,
                            chrom in enumerate(fd.references)}

    if bam:
        with pysam.AlignmentFile(bam, "rb") as fd:
            chroms_order = {chrom: chrom_index for chrom_index,
                            chrom in enumerate(fd.references)}

    return chroms_order


def safe_read_info_dict(d, field, t=str, default_val=""):
    return t(d[field]) if field in d else default_val


def prob2phred(p, max_phred=100):
    '''Convert prob to Phred-scale quality score.'''
    assert 0 <= p <= 1
    if p == 1:
        Q = max_phred
    elif p == 0:
        Q = 0
    else:
        Q = -10 * np.log10(1 - p)
        if Q > max_phred:
            Q = max_phred
    return Q
