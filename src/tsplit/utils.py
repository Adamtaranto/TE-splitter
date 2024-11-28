from datetime import datetime, timezone
import logging
import os
import re
import shutil
import sys

from Bio import SeqIO


def check_tools(required_tools=[], optional_tools=[]):
    """
    Check if required and optional tools are available on the system's PATH.

    Args:
        required_tools (list): List of required tool names.
        optional_tools (list): List of optional tool names.

    Raises:
        RuntimeError: If any required tool is not found.
    """
    missing_required_tools = []

    def print_message(tool, path, color_code):
        """
        Print a message to stderr with the tool name and path in the specified color.

        Args:
            tool (str): The name of the tool.
            path (str): The path to the tool.
            color_code (str): The ANSI color code for the message.
        """
        tool_padded = tool.ljust(15)
        if path:
            message = f"{color_code}{tool_padded}\t{path}\033[0m"
        else:
            message = f"{color_code}{tool_padded}\tNOT FOUND\033[0m"
        print(message, file=sys.stderr)

    # Check required tools
    logging.info("Checking for dependencies:")
    for tool in required_tools:
        path = shutil.which(tool)
        if path:
            print_message(tool, path, "\033[92m")  # Green
        else:
            print_message(tool, None, "\033[91m")  # Red
            missing_required_tools.append(tool)

    # Check optional tools
    for tool in optional_tools:
        path = shutil.which(tool)
        if path:
            print_message(tool, path, "\033[92m")  # Green
        else:
            print_message(tool, None, "\033[93m")  # Yellow

    # Raise error if any required tool is missing
    if missing_required_tools:
        error_message = "ERROR: Some required tools could not be found: " + ", ".join(
            missing_required_tools
        )
        logging.error(error_message)
        raise RuntimeError(
            "Missing required tools: " + ", ".join(missing_required_tools)
        )


def tSplitchecks(args):
    """
    Perform housekeeping tasks: Create output files/dirs and temp dirs as required.

    Args:
        args: Argument parser object containing input parameters.

    Returns:
        str: Full path to the output file.

    Raises:
        FileNotFoundError: If the input file does not exist.
    """
    # Check if the input file exists
    if not os.path.isfile(args.infile):
        logging.error("Input sequence file does not exist. Quitting.")
        raise FileNotFoundError(f"Input sequence file '{args.infile}' does not exist.")
    logging.info(f"Input file found: {args.infile}")

    # Determine the output directory
    if args.outdir:
        absOutDir = os.path.abspath(args.outdir)
        if not os.path.isdir(absOutDir):
            os.makedirs(absOutDir)
            logging.info(f"Creating output directory: {absOutDir}")
        outDir = absOutDir
    else:
        outDir = os.getcwd()
    logging.debug(f"Set output directory: {outDir}")

    # Set the prefix for output files
    if not args.prefix:
        prefix = os.path.splitext(os.path.basename(args.infile))[0]
    else:
        prefix = args.prefix
    logging.debug(f"Set prefix: {prefix}")

    # Create the output file path
    outfile = prefix + "_tsplit_output.fasta"
    outpath = os.path.join(outDir, outfile)
    logging.debug(f"Set outfile target: {outpath}")

    # Return the full path to the output file
    return outpath


def cleanID(s):
    """
    Remove non alphanumeric characters from string.
    Replace whitespace with underscores.
    """
    s = re.sub(r"[^\w\s]", "", s)
    s = re.sub(r"\s+", "_", s)
    return s


def segWrite(outfile, segs=None):
    """
    Take a generator object yielding seqrecords and
    write each to outfile in fasta format.
    """
    seqcount = 0
    if segs:
        with open(outfile, "w") as handle:
            for seq in segs:
                seqcount += 1
                SeqIO.write(seq, handle, "fasta")
        if seqcount == 0:
            os.remove(outfile)
