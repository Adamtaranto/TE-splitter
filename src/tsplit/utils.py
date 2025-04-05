"""
Utility functions for tSplit.

This module provides helper functions for file operations, dependency checking,
sequence ID handling, and other utility tasks used across the TE-splitter application.
"""

import logging
import os
import re
import shutil
import sys
from typing import Any, Generator, List, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def check_tools(
    required_tools: Optional[List[str]] = None,
    optional_tools: Optional[List[str]] = None,
) -> None:
    """
    Check if required and optional tools are available on the system's PATH.

    Parameters
    ----------
    required_tools : list of str, optional
        List of required tool names.
    optional_tools : list of str, optional
        List of optional tool names.

    Raises
    ------
    RuntimeError
        If any required tool is not found.

    Notes
    -----
    Tools in the required_tools list will cause the program to exit if not found.
    Tools in the optional_tools list will only produce a warning if not found.
    """
    # Initialize empty lists if None is provided
    if required_tools is None:
        required_tools = []
    if optional_tools is None:
        optional_tools = []

    missing_required_tools = []

    def print_message(tool: str, path: Optional[str], color_code: str) -> None:
        """
        Print a message to stderr with the tool name and path in the specified color.

        Parameters
        ----------
        tool : str
            The name of the tool.
        path : str or None
            The path to the tool, or None if not found.
        color_code : str
            The ANSI color code for the message.
        """
        tool_padded = tool.ljust(15)
        if path:
            message = f'{color_code}{tool_padded}\t{path}\033[0m'
        else:
            message = f'{color_code}{tool_padded}\tNOT FOUND\033[0m'
        print(message, file=sys.stderr)

    # Check required tools
    logging.info('Checking for dependencies:')
    for tool in required_tools:
        path = shutil.which(tool)  # Get the full path of the executable
        if path:
            print_message(tool, path, '\033[92m')  # Green for found tools
        else:
            print_message(tool, None, '\033[91m')  # Red for missing required tools
            missing_required_tools.append(tool)

    # Check optional tools
    for tool in optional_tools:
        path = shutil.which(tool)
        if path:
            print_message(tool, path, '\033[92m')  # Green for found tools
        else:
            print_message(tool, None, '\033[93m')  # Yellow for missing optional tools

    # Raise error if any required tool is missing
    if missing_required_tools:
        error_message = 'ERROR: Some required tools could not be found: ' + ', '.join(
            missing_required_tools
        )
        logging.error(error_message)
        raise RuntimeError(
            'Missing required tools: ' + ', '.join(missing_required_tools)
        )


def tSplitchecks(args: Any) -> str:
    """
    Perform housekeeping tasks for TE-splitter.

    Creates output directories, sets file naming conventions,
    and validates input files.

    Parameters
    ----------
    args : argparse.Namespace
        Argument parser object containing input parameters.

    Returns
    -------
    str
        Full path to the output file.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    """
    # Check if the input file exists
    if not os.path.isfile(args.infile):
        logging.error('Input sequence file does not exist. Quitting.')
        raise FileNotFoundError(f"Input sequence file '{args.infile}' does not exist.")
    logging.info(f'Input file found: {args.infile}')

    # Determine the output directory - create if it doesn't exist
    if args.outdir:
        absOutDir = os.path.abspath(args.outdir)
        if not os.path.isdir(absOutDir):
            os.makedirs(absOutDir)
            logging.info(f'Creating output directory: {absOutDir}')
        outDir = absOutDir
    else:
        outDir = os.getcwd()  # Use current working directory if none specified
    logging.debug(f'Set output directory: {outDir}')

    # Set the prefix for output files - use input filename if none provided
    if not args.prefix:
        prefix = os.path.splitext(os.path.basename(args.infile))[0]
    else:
        prefix = args.prefix
    logging.debug(f'Set prefix: {prefix}')

    # Create the output file path by combining directory and filename
    outfile = prefix + '_tsplit_output.fasta'
    outpath = os.path.join(outDir, outfile)
    logging.debug(f'Set outfile target: {outpath}')

    # Return the full path to the output file
    return outpath


def cleanID(s: str) -> str:
    """
    Remove non-alphanumeric characters from string and replace whitespace.

    Parameters
    ----------
    s : str
        Input string to be cleaned.

    Returns
    -------
    str
        Cleaned string with special characters removed and
        whitespace replaced with underscores.

    Notes
    -----
    This function is useful for ensuring sequence IDs are compatible
    with various bioinformatics tools.
    """
    # Remove any character that isn't alphanumeric or whitespace
    s = re.sub(r'[^\w\s]', '', s)
    # Replace any whitespace sequence with a single underscore
    s = re.sub(r'\s+', '_', s)
    return s


def segWrite(
    outfile: str, segs: Optional[Generator[SeqRecord, None, None]] = None
) -> None:
    """
    Write sequence records to a FASTA file.

    Takes a generator object yielding SeqRecord objects and
    writes each to the specified output file in FASTA format.
    If no sequences are written, the output file is removed.

    Parameters
    ----------
    outfile : str
        Path to the output file.
    segs : generator of Bio.SeqRecord.SeqRecord, optional
        Generator yielding SeqRecord objects to write.

    Returns
    -------
    None
        This function writes data to a file and doesn't return any value.

    Notes
    -----
    If no sequences are written (empty generator), the output file
    will be automatically removed to avoid empty files.
    """
    seqcount = 0
    if segs:
        # Open the output file for writing
        with open(outfile, 'w') as handle:
            for seq in segs:
                seqcount += 1
                SeqIO.write(seq, handle, 'fasta')

        # Clean up empty files to avoid clutter
        if seqcount == 0:
            os.remove(outfile)
