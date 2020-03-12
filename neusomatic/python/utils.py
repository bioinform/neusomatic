#-------------------------------------------------------------------------
# utils.py
# Utility functions
#-------------------------------------------------------------------------
import os
import shutil
import shlex
import subprocess
import logging
import traceback
import tempfile

import pysam
import numpy as np


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


def run_shell_command(command, stdout=None, stderr=None, run_logger=None):
    stdout_fd = open(stdout, "w") if stdout else None
    stderr_fd = open(stderr, "w") if stderr else None
    my_logger = logger
    if run_logger:
        my_logger = run_logger

    fixed_command = shlex.split(command)
    my_logger.info("Running command: {}".format(fixed_command))
    returncode = subprocess.check_call(
        fixed_command, stdout=stdout_fd, stderr=stderr_fd)
    if stdout_fd:
        stdout_fd.close()
    if stderr_fd:
        stderr_fd.close()
    return returncode


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


def run_bedtools_cmd(command, output_fn=None, run_logger=None):
    if output_fn is None:
        tmpfn = tempfile.NamedTemporaryFile(
            prefix="tmpbed_", suffix=".bed", delete=False)
        output_fn = tmpfn.name
    stderr_file = output_fn + ".stderr"
    if run_logger is None:
        run_logger = logger
    try:
        returncode = run_shell_command(command, stdout=output_fn, stderr=stderr_file,
                                       run_logger=run_logger)
        os.remove(stderr_file)
        return output_fn
    except Exception as ex:
        run_logger.error(traceback.format_exc())
        run_logger.error(ex)
        err_msg = "Command {} failed.".format(command)
        if os.path.exists(stderr_file) and os.path.getsize(stderr_file):
            err_msg += "\nPlease check error log at {}".format(stderr_file)
        raise Exception(err_msg)


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
