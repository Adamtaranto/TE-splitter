from shlex import quote
import logging
import subprocess
import sys


def run_cmd(cmd, verbose=False, workingDir=None):
    """
    Execute a command in the specified directory.

    Args:
        cmd (str): The command to execute.
        verbose (bool): Whether to print verbose output.
        workingDir (str): The directory in which to execute the command.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit code.
    """
    if verbose:
        print("\nRunning command:", cmd, flush=True)
    try:
        output = subprocess.check_output(
            cmd, shell=True, stderr=subprocess.STDOUT, cwd=workingDir
        )
        if verbose:
            print(output.decode())
    except subprocess.CalledProcessError as error:
        print(
            "The following command failed with exit code",
            error.returncode,
            file=sys.stderr,
        )
        print(cmd, file=sys.stderr)
        print("\nThe output was:\n", file=sys.stderr)
        print(error.output.decode(), file=sys.stderr)
        logging.error(f"Error running command: {cmd}")
        raise


def makeBlast(seq=None, outfile=None, pid=60):
    """
    Construct the blastn command.
    """
    cmd = (
        'blastn -word_size 4 -outfmt "6 qstart qend sstart send length positive pident qlen slen qframe sframe qseqid sseqid" -query '
        + quote(str(seq))
        + " -subject "
        + quote(str(seq))
        + " -out "
        + quote(str(outfile))
        + " -perc_identity "
        + str(pid)
    )
    return cmd


"""
Recreate pymummer-like coords file output from blast.

#Blast field, coords field, Description

qstart			[S1] 	Start of the alignment region in the reference sequence 
qend			[E1] 	End of the alignment region in the reference sequence 
sstart			[S2] 	Start of the alignment region in the query sequence 
send 			[E2] 	End of the alignment region in the query sequence 
length			[LEN 1] Length of the alignment region in the reference sequence
positive		[LEN 2] Length of the alignment region in the query sequence (#"positive" is just a filler as blast won't all repeated fields)
pident			[% IDY] Percent identity of the alignment 
qlen			[LEN R] Length of the reference sequence 
slen			[LEN Q] Length of the query sequence 
qframe sframe	[FRM] 	Reading frame for the reference AND query sequence alignments respectively 
qseqid sseqid	[TAGS] 	The reference AND query FastA IDs respectively. All output coordinates and lengths are relative to the forward strand of the reference DNA sequence.
"""
