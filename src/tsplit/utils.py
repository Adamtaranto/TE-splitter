from datetime import datetime
import os
import sys
from Bio import SeqIO
from collections import Counter
import re
import shutil
import logging


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
        if path:
            message = f"{color_code}Tool found: {tool} at {path}\033[0m"
        else:
            message = f"{color_code}Tool not found: {tool}\033[0m"
        print(message, file=sys.stderr)

    # Check required tools
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
        error_message = (
            "ERROR: Some required tools could not be found: "
            + ", ".join(missing_required_tools)
        )
        logging.error(error_message)
        raise RuntimeError("Missing required tools: " + ", ".join(missing_required_tools))



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
            logging.info(f"Output directory created: {absOutDir}")
        outDir = absOutDir
    else:
        outDir = os.getcwd()
    logging.info(f"Output directory set to: {outDir}")

    # Set the prefix for output files
    if not args.prefix:
        prefix = os.path.splitext(os.path.basename(args.infile))[0]
    else:
        prefix = args.prefix
    logging.info(f"Prefix set to: {prefix}")

    # Create the output file path
    outfile = prefix + "_tsplit_output.fasta"
    outpath = os.path.join(outDir, outfile)
    logging.info(f"Output file name set to: {outpath}")

    # Return the full path to the output file
    return outpath


def getTimestring():
    """
    Return int only string of current datetime with milliseconds.
    """
    (dt, micro) = datetime.utcnow().strftime("%Y%m%d%H%M%S.%f").split(".")
    dt = "%s%03d" % (dt, int(micro) / 1000)
    return dt



def cleanID(s):
    """
    Remove non alphanumeric characters from string.
    Replace whitespace with underscores.
    """
    s = re.sub(r"[^\w\s]", "", s)
    s = re.sub(r"\s+", "_", s)
    return s

## Fix: Do not load fasta into memory!
def importFasta2List(file):
    """Load elements from multifasta file. Check that seq IDs are unique."""
    # Read in elements from multifasta file, convert seqrecord iterator to list
    records = list(SeqIO.parse(file, "fasta"))
    # Check names are unique
    checkUniqueID(records)
    # If unique, return record list.
    return records



# Fix: Do not load fasta into genome!
def checkUniqueID(records):
    """
    Check that IDs for input elements are unique.
    """
    seqIDs = [records[x].id for x in range(len(records))]
    IDcounts = Counter(seqIDs)
    duplicates = [k for k, v in IDcounts.items() if v > 1]
    if duplicates:
        print("Input sequence IDs not unique. Quiting.")
        print(duplicates)
        sys.exit(1)
    else:
        pass

## Fix: Do not load fasta into memory!
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
            
# Fix: Do not load fasta into genome!
def importFasta(file):
    """
    Load elements from multifasta file. Check that seq IDs are unique.
    """
    # Read in elements from multifasta file, convert seqrecord iterator to list
    records = list(SeqIO.parse(file, "fasta"))
    # Check names are unique
    checkUniqueID(records)
    # If unique, return records as dict keyed by seq id
    recordsDict = dict()
    for rec in records:
        recordsDict[rec.id] = rec
    return recordsDict
