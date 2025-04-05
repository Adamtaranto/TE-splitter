"""
Wrapper functions for external tool execution.

This module provides utility functions for executing external commands
and constructing specific command strings for tools like BLASTN.
"""

import logging
from shlex import quote
import subprocess
import sys
from typing import Optional


def run_cmd(
    cmd: str, verbose: bool = False, workingDir: Optional[str] = None
) -> Optional[str]:
    """
    Execute a command in the specified directory.

    Parameters
    ----------
    cmd : str
        The command to execute.
    verbose : bool, optional
        Whether to print verbose output. Default is False.
    workingDir : str, optional
        The directory in which to execute the command. Default is None.

    Returns
    -------
    str
        The decoded output of the command if successful.

    Raises
    ------
    subprocess.CalledProcessError
        If the command returns a non-zero exit code.
    """
    if verbose:
        print('\nRunning command:', cmd, flush=True)
    try:
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT, cwd=workingDir
        )
        if verbose:
            print(output.decode())
        return output.decode()
    except subprocess.CalledProcessError as error:
        print(
            'The following command failed with exit code',
            error.returncode,
            file=sys.stderr,
        )
        print(cmd, file=sys.stderr)
        print('\nThe output was:\n', file=sys.stderr)
        print(error.output.decode(), file=sys.stderr)
        logging.error(f'Error running command: {cmd}')
        raise


def makeBlast(
    seq: Optional[str] = None, outfile: Optional[str] = None, pid: int = 60
) -> str:
    """
    Construct the blastn command.

    Creates a command string to run BLASTN for self-alignment with specific
    formatting options.

    Parameters
    ----------
    seq : str, optional
        Path to the sequence file to be used as both query and subject.
    outfile : str, optional
        Path where the BLAST results will be saved.
    pid : int, optional
        Minimum percent identity threshold. Default is 60.

    Returns
    -------
    str
        The constructed blastn command string.

    Notes
    -----
    The BLAST output format is tab-delimited with the following fields:
    qstart, qend, sstart, send, length, positive, pident, qlen, slen, qframe,
    sframe, qseqid, sseqid

    These fields correspond to coordinates file output fields:
    - qstart [S1]: Start of alignment region in reference sequence
    - qend [E1]: End of alignment region in reference sequence
    - sstart [S2]: Start of alignment region in query sequence
    - send [E2]: End of alignment region in query sequence
    - length [LEN 1]: Length of alignment region in reference sequence
    - positive [LEN 2]: Length of alignment region in query sequence
    - pident [% IDY]: Percent identity of the alignment
    - qlen [LEN R]: Length of reference sequence
    - slen [LEN Q]: Length of query sequence
    - qframe sframe [FRM]: Reading frames for reference and query
    - qseqid sseqid [TAGS]: Reference and query FastA IDs
    """
    cmd = (
        'blastn -word_size 4 -outfmt "6 qstart qend sstart send length positive pident qlen slen qframe sframe qseqid sseqid" -query '
        + quote(str(seq))
        + ' -subject '
        + quote(str(seq))
        + ' -out '
        + quote(str(outfile))
        + ' -perc_identity '
        + str(pid)
    )
    return cmd
