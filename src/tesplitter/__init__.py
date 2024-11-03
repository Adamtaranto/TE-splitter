from collections import Counter
from datetime import datetime
import argparse
import os
import re
import shutil
import sys

from Bio import SeqIO
from pymummer import coords_file, alignment, nucmer

from .runBlastn import makeBlast, run_blast

def dochecks(args):
    """Housekeeping tasks: Create output files/dirs and temp dirs as required."""
    if not os.path.isfile(args.infile):
        print("Input sequence file does not exist. Quitting.")
        sys.exit(1)
    # Make outDir if does not exist else set to current dir.
    if args.outdir:
        absOutDir = os.path.abspath(args.outdir)
        if not os.path.isdir(args.outdir):
            os.makedirs(absOutDir)
        outDir = absOutDir
    else:
        outDir = os.getcwd() 
    # Make temp directory
    tempDir = os.path.join(os.getcwd(),"temp_" + getTimestring())
    os.makedirs(tempDir)
    # Set prefix to infile basename if none
    if not args.prefix:
        prefix = os.path.splitext(os.path.basename(args.infile))[0]
    else:
        prefix = args.prefix
    # Create outfile paths
    outfile = prefix + "_tsplit_output.fasta"
    outpath = os.path.join(outDir,outfile)
    # Return full path to output file and temp directory
    return outpath,tempDir

def getTimestring():
    """Return int only string of current datetime with milliseconds."""
    (dt, micro) = datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.')
    dt = "%s%03d" % (dt, int(micro) / 1000)
    return dt

def checkUniqueID(records):
    """Check that IDs for input elements are unique."""
    seqIDs = [records[x].id for x in range(len(records))]
    IDcounts = Counter(seqIDs)
    duplicates = [k for k, v in IDcounts.items() if v > 1]
    if duplicates:
        print("Input sequence IDs not unique. Quiting.")
        print(duplicates)
        sys.exit(1)
    else:
        pass

def manageTemp(record=None, tempPath=None, scrub=False):
    """Create single sequence fasta files or scrub temp files."""
    if scrub and tempPath:
        try:
            os.remove(tempPath)
        except OSError:
            pass
    else:
        with open(tempPath, "w") as f:
            SeqIO.write(record, f, "fasta")

def importFasta(file):
    """Load elements from multifasta file. Check that seq IDs are unique."""
    # Read in elements from multifasta file, convert seqrecord iterator to list.
    records = list(SeqIO.parse(file, "fasta"))
    # Check names are unique
    checkUniqueID(records)
    # If unique, return record list.
    return records

def cleanID(s):
    """Remove non alphanumeric characters from string. Replace whitespace with underscores."""
    s = re.sub(r"[^\w\s]", '', s)
    s = re.sub(r"\s+", '_', s)
    return s

def getLTRs(elements=None, flankdist=10, minid=80, minterm=10, minseed=5, diagfactor=0.3, report='split', temp=None,keeptemp=False,alignTool='nucmer',verbose=False):
    """Align elements to self and attempt to identify LTRs."""
    # Set temp directory to cwd if none.
    if not temp:
        temp = os.getcwd()
    # For each candidate LTR element
    for rec in elements:
        # Create temp paths for single element fasta and alignment coords
        tempFasta = os.path.join(temp, cleanID(rec.id) + '.fasta')
        tempCoords = tempFasta + '.coords'
        # Write current element to single fasta
        manageTemp(record=rec, tempPath=tempFasta, scrub=False)
        # Align to self with nucmer
        if alignTool == 'nucmer':
            # Compose Nucmer script for current element vs self
            runner = nucmer.Runner(	tempFasta, tempFasta, tempCoords,
                                    min_id		=	minid, 
                                    min_length	=	minseed,
                                    diagfactor	=	diagfactor,
                                    mincluster	=	minterm,
                                    breaklen	=	200,
                                    maxmatch	=	True,
                                    simplify	=	False
                                    )
            # Execute nucmer
            runner.run()
        elif alignTool == 'blastn':
            # Alternatively, use blastn as search tool and write nucmer.coords-like output.
            cmds = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
            run_blast(cmds, verbose=verbose)
        # Import coords file to iterator object
        file_reader = coords_file.reader(tempCoords)
        # Exclude hits to self. Also converts iterator output to stable list
        alignments = [hit for hit in file_reader if not hit.is_self_hit()]
        # Filter hits less than min length (Done internally for nucmer, not blastn.)
        alignments = [hit for hit in alignments if hit.ref_end - hit.ref_start >= minterm]
        # Filter for hits on same strand i.e. tandem repeats / LTRs
        alignments = [hit for hit in alignments if hit.on_same_strand()]
        # Filter for 5' repeats which begin within x bases of element start
        alignments = [hit for hit in alignments if hit.ref_start <= flankdist]
        # Filter for 5' repeats whose 3' match ends within x bases of element end
        alignments = [hit for hit in alignments if len(rec) - hit.qry_end <= flankdist]
        # Scrub overlappying ref / query segments, and also complementary 3' to 5' flank hits
        alignments = [hit for hit in alignments if hit.ref_end < hit.qry_start]
        # Sort largest to smallest dist between end of ref (subject) and start of query (hit)
        # x.qry_start (3') - x.ref_end (5') = Length of internal segment
        alignments = sorted(alignments, key=lambda x: (x.qry_start - x.ref_end), reverse=True)
        # If alignments exist after filtering report features using alignment pair with largest 
        # internal segment i.e. first element in sorted list.
        if alignments:
            if verbose:
                [print(x) for x in alignments]
            if report == 'all':
                # yield original element
                yield rec
            if report in ['split','external']:
                # yield LTR slice - append "_LTR"
                extSeg = rec[alignments[0].ref_start:alignments[0].ref_end + 1]
                extSeg.id = extSeg.id + "_LTR"
                extSeg.name = extSeg.id
                extSeg.description = "[" + rec.id + " LTR segment]"
                yield extSeg
            if report in ['split','internal']:
                # yield internal slice - append "_I"
                intSeg = rec[alignments[0].ref_end:alignments[0].qry_start + 1] 
                intSeg.id = intSeg.id + "_I"
                intSeg.name = intSeg.id
                intSeg.description = "[" + rec.id + " internal segment]"
                yield intSeg
        else:
            # If alignment list empty after filtering print alert and continue
            print('No LTRs found for candidate element: %s' % rec.id)
        # Scrub single fasta and coords file for current element.
        if not keeptemp:
            manageTemp(tempPath=tempFasta, scrub=True)
            manageTemp(tempPath=tempCoords, scrub=True)

def getTIRs(elements=None, flankdist=10, minid=80, minterm=10, minseed=5, diagfactor=0.3, mites=False, report='split', temp=None,keeptemp=False,alignTool='nucmer',verbose=False):
    """ 
    Align elements to self and attempt to identify TIRs. 
    Optionally attempt to construct synthetic MITEs from TIRs.
    """
    # Set temp directory to cwd if none.
    if not temp:
        temp = os.getcwd()
    # For each candidate LTR element
    for rec in elements:
        # Create temp paths for single element fasta and alignment coords
        tempFasta = os.path.join(temp, cleanID(rec.id) + '.fasta')
        tempCoords = tempFasta + '.coords'
        # Write current element to single fasta
        manageTemp(record=rec, tempPath=tempFasta, scrub=False)
        # Align to self with nucmer
        if alignTool == 'nucmer':
            # Compose Nucmer script for current element vs self
            runner = nucmer.Runner(	tempFasta, tempFasta, tempCoords,
                                    min_id		=	minid,
                                    min_length	=	minseed,
                                    diagfactor	=	diagfactor,
                                    mincluster	=	minterm,
                                    breaklen	=	200,
                                    maxmatch	=	True,
                                    simplify	=	False
                                    )
            # Execute nucmer
            runner.run()
        elif alignTool == 'blastn':
            # Alternatively, use blastn as search tool and write nucmer.coords-like output.
            cmds = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
            run_blast(cmds,verbose=verbose)
        # Import coords file to iterator object
        file_reader = coords_file.reader(tempCoords)
        # Exclude hits to self. Also converts iterator output to stable list
        alignments = [hit for hit in file_reader if not hit.is_self_hit()]
        # Filter hits less than min length (Done internally for nucmer, not blastn.)
        alignments = [hit for hit in alignments if hit.ref_end - hit.ref_start >= minterm]
        # Filter for hits on same strand i.e. tandem repeats / LTRs
        alignments = [hit for hit in alignments if not hit.on_same_strand()]
        # Filter for 5' repeats which begin within x bases of element start
        alignments = [hit for hit in alignments if hit.ref_start <= flankdist]
        # Scrub overlappying ref / query segments, and also complementary 3' to 5' flank hits
        alignments = [hit for hit in alignments if hit.ref_end < hit.qry_end]
        # Sort largest to smallest dist between end of ref (subject) and start of query (hit)
        # x.qry_start - x.ref_end = length of internal segment
        alignments = sorted(alignments, key=lambda x: (x.qry_end - x.ref_end), reverse=True)
        # If alignments exist after filtering report features using alignment pair with largest 
        # internal segment i.e. first element in sorted list.
        if alignments:
            if verbose:
                [print(x) for x in alignments]
            if report == 'all':
                yield rec
            if report in ['split','external']:
                # yield TIR slice - append "_TIR"
                extSeg = rec[alignments[0].ref_start:alignments[0].ref_end + 1]
                extSeg.id = extSeg.id + "_TIR"
                extSeg.name = extSeg.id
                extSeg.description = "[" + rec.id + " TIR segment]"
                yield extSeg
            if report in ['split','internal']:
                # yield internal slice - append "_I"
                intSeg = rec[alignments[0].ref_end:alignments[0].qry_end + 1]
                intSeg.id = intSeg.id + "_I"
                intSeg.name = intSeg.id
                intSeg.description = "[" + rec.id + " internal segment]"
                yield intSeg
            if mites:
                # Assemble TIRs into hypothetical MITEs
                synMITE = rec[alignments[0].ref_start:alignments[0].ref_end + 1] + rec[alignments[0].qry_end:alignments[0].qry_start + 1]
                synMITE.id = synMITE.id + "_synMITE"
                synMITE.name = synMITE.id
                synMITE.description = "[Synthetic MITE constructed from " + rec.id + " TIRs]"
                yield synMITE
        else:
            # If alignment list empty after filtering print alert and continue
            print('No TIRs found for candidate element: %s' % rec.id)
        # Scrub single fasta and coords file for current element.
        if not keeptemp:
            manageTemp(tempPath=tempFasta, scrub=True)
            manageTemp(tempPath=tempCoords, scrub=True)

def segWrite(outfile,segs=None):
    """Take a generator object yielding seqrecords and write each to outfile in fasta format."""
    seqcount = 0
    if segs:
        with open(outfile, "w") as handle:
            for seq in segs:
                seqcount += 1
                SeqIO.write(seq, handle, "fasta")
        if seqcount == 0:
            os.remove(outfile)

"""
# Useful attributes of pymummer objects:
[x.ref_start for x in alignments]
[x.ref_end for x in alignments]
[x.qry_start for x in alignments]
[x.qry_end for x in alignments]
[x.hit_length_ref for x in alignments]
[x.hit_length_qry for x in alignments]
[x.percent_identity for x in alignments]
[x.ref_length for x in alignments]
[x.qry_length for x in alignments]
[x.frame for x in alignments]
[x.ref_name for x in alignments]
[x.qry_name for x in alignments]

## Don't use these, bizzaro format. Not indexed to 0. Cannot sort as ints.
#coord.reverse_query()
#coord.reverse_reference()
#coord.qry_coords()
#coord.ref_coords()
"""