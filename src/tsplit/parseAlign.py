import os
import logging
import tempfile
import shutil
from pymummer import coords_file, nucmer
from Bio import SeqIO
from tsplit.wrapping import run_cmd, makeBlast
from tsplit.utils import cleanID, getTimestring

def getTIRs(
    elements=[],
    flankdist=2,
    minid=80,
    minterm=10,
    minseed=5,
    diagfactor=0.3,
    mites=False,
    report="split",
    temp=None,
    keeptemp=False,
    alignTool="nucmer",
    verbose=True,
):
    """
    Align elements to self and attempt to identify TIRs.
    Optionally attempt to construct synthetic MITEs from TIRs.

    Args:
        elements (list): List of sequence records to be processed.
        flankdist (int): Maximum distance from element start for TIR candidates.
        minid (float): Minimum identity between terminal repeat pairs.
        minterm (int): Minimum length for a terminal repeat to be considered.
        minseed (int): Minimum seed length for nucmer.
        diagfactor (float): Diagonal factor for nucmer.
        mites (bool): Whether to attempt to construct synthetic MITEs.
        report (str): Reporting mode for TIRs.
        temp (str): Path to the temporary directory.
        keeptemp (bool): Whether to keep the temporary directory after processing.
        alignTool (str): Alignment tool to use ('nucmer' or 'blastn').
        verbose (bool): Whether to print verbose output.

    Yields:
        SeqRecord: Segments of the sequence based on the reporting mode.
    """
    # Set temp directory to cwd if none is provided
    if not temp:
        temp = os.getcwd()
    
    tempName = os.path.join(temp, f"tsplit_temp_{getTimestring()}")
    
    # Create a unique temporary directory
    with tempfile.TemporaryDirectory(dir=tempName) as tempDir:
        logging.info(f"Temporary directory created: {tempDir}")

        try:
            # For each candidate TIR element
            for rec in elements:
                # Create temp paths for single element fasta and alignment coords
                tempFasta = os.path.join(tempDir, cleanID(rec.id) + ".fasta")
                tempCoords = os.path.join(tempDir, cleanID(rec.id) + "_" + alignTool + ".coords")
                # Write current element to single fasta
                with open(tempFasta, "w") as f:
                    SeqIO.write(rec, f, "fasta")
                # Align to self with nucmer
                if alignTool == "nucmer":
                    # Compose Nucmer script for current element vs self
                    runner = nucmer.Runner(
                        tempFasta,
                        tempFasta,
                        tempCoords,
                        min_id=minid,
                        min_length=minseed,
                        diagfactor=diagfactor,
                        mincluster=minterm,
                        breaklen=200,
                        maxmatch=True,
                        simplify=False,
                    )
                    # Execute nucmer
                    runner.run()
                elif alignTool == "blastn":
                    # Alternatively, use blastn as search tool and write nucmer.coords-like output.
                    cmd = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
                    run_cmd(cmd, verbose=verbose, tempDir=tempDir)
                # Import coords file to iterator object
                file_reader = coords_file.reader(tempCoords)
                # Exclude hits to self. Also converts iterator output to stable list
                alignments = [hit for hit in file_reader if not hit.is_self_hit()]
                # Filter hits less than min length (Done internally for nucmer, not blastn.)
                alignments = [
                    hit for hit in alignments if hit.ref_end - hit.ref_start >= minterm
                ]
                # Filter for hits on same strand i.e. tandem repeats / LTRs
                alignments = [hit for hit in alignments if not hit.on_same_strand()]
                # Filter for 5' repeats which begin within x bases of element start
                alignments = [hit for hit in alignments if hit.ref_start <= flankdist]
                # Scrub overlapping ref / query segments, and also complementary
                # 3' to 5' flank hits
                alignments = [hit for hit in alignments if hit.ref_end < hit.qry_end]
                # Sort largest to smallest dist between end of ref (subject) and start
                # of query (hit)
                # x.qry_end - x.ref_end = 5'end of right TIR - 3' end of left
                # TIR = length of internal segment
                # TIR pair with smallest internal segment (longest TIRs) is first in list.
                alignments = sorted(
                    alignments, key=lambda x: (x.qry_end - x.ref_end), reverse=False
                )
                # If alignments exist after filtering report features using alignment
                # pair with largest internal segment i.e. first element in sorted list.
                if alignments:
                    if verbose:
                        [print(x) for x in alignments]
                    if report in ["split", "external", "all"]:
                        # yield TIR slice - append "_TIR"
                        extSeg = rec[alignments[0].ref_start : alignments[0].ref_end + 1]
                        extSeg.id = extSeg.id + "_TIR"
                        extSeg.name = extSeg.id
                        extSeg.description = "[" + rec.id + " TIR segment]"
                        yield extSeg
                    if report in ["split", "internal", "all"]:
                        # yield internal slice - append "_I"
                        intSeg = rec[alignments[0].ref_end : alignments[0].qry_end + 1]
                        intSeg.id = intSeg.id + "_I"
                        intSeg.name = intSeg.id
                        intSeg.description = "[" + rec.id + " internal segment]"
                        yield intSeg
                    if report == "all":
                        yield rec
                    if mites:
                        # Assemble TIRs into hypothetical MITEs
                        synMITE = (
                            rec[alignments[0].ref_start : alignments[0].ref_end + 1]
                            + rec[alignments[0].qry_end : alignments[0].qry_start + 1]
                        )
                        synMITE.id = synMITE.id + "_synMITE"
                        synMITE.name = synMITE.id
                        synMITE.description = (
                            "[Synthetic MITE constructed from " + rec.id + " TIRs]"
                        )
                        yield synMITE
                else:
                    # If alignment list empty after filtering, print alert and continue
                    logging.info(f"No TIRs found for candidate element: {rec.id}")
        finally:
            # Clean up the temporary directory if keeptemp is False
            if keeptemp:
                shutil.copytree(tempDir, os.path.join(temp, "kept_temp"))
                logging.info(f"Temporary directory retained: {tempDir}")
            else:
                logging.info(f"Temporary directory deleted: {tempDir}")

#def getLTRs(elements=None, flankdist=10, minid=80, minterm=10, minseed=5, diagfactor=0.3, report='split', temp=None,keeptemp=False,alignTool='nucmer',verbose=True):
#    """Align elements to self and attempt to identify LTRs."""
#    # Set temp directory to cwd if none.
#    if not temp:
#        temp = os.getcwd()
#    # For each candidate LTR element
#    for rec in elements:
#        # Create temp paths for single element fasta and alignment coords
#        tempFasta = os.path.join(temp, cleanID(rec.id) + '.fasta')
#        tempCoords = tempFasta + '.coords'
#        # Write current element to single fasta
#        manageTemp(record=rec, tempPath=tempFasta, scrub=False)
#        # Align to self with nucmer
#        if alignTool == 'nucmer':
#            # Compose Nucmer script for current element vs self
#            runner = nucmer.Runner(	tempFasta, tempFasta, tempCoords,
#                                    min_id		=	minid, 
#                                    min_length	=	minseed,
#                                    diagfactor	=	diagfactor,
#                                    mincluster	=	minterm,
#                                    breaklen	=	200,
#                                    maxmatch	=	True,
#                                    simplify	=	False
#                                    )
#            # Execute nucmer
#            runner.run()
#        elif alignTool == 'blastn':
#            # Alternatively, use blastn as search tool and write nucmer.coords-like output.
#            cmds = makeBlast(seq=tempFasta, outfile=tempCoords, pid=minid)
#            run_cmd(cmds, verbose=verbose)
#        # Import coords file to iterator object
#        file_reader = coords_file.reader(tempCoords)
#        # Exclude hits to self. Also converts iterator output to stable list
#        alignments = [hit for hit in file_reader if not hit.is_self_hit()]
#        # Filter hits less than min length (Done internally for nucmer, not blastn.)
#        alignments = [hit for hit in alignments if hit.ref_end - hit.ref_start >= minterm]
#        # Filter for hits on same strand i.e. tandem repeats / LTRs
#        alignments = [hit for hit in alignments if hit.on_same_strand()]
#        # Filter for 5' repeats which begin within x bases of element start
#        alignments = [hit for hit in alignments if hit.ref_start <= flankdist]
#        # Filter for 5' repeats whose 3' match ends within x bases of element end
#        alignments = [hit for hit in alignments if len(rec) - hit.qry_end <= flankdist]
#        # Scrub overlappying ref / query segments, and also complementary 3' to 5' flank hits
#        alignments = [hit for hit in alignments if hit.ref_end < hit.qry_start]
#        # Sort largest to smallest dist between end of ref (subject) and start of query (hit)
#        # x.qry_start (3') - x.ref_end (5') = Length of internal segment
#        alignments = sorted(alignments, key=lambda x: (x.qry_start - x.ref_end), reverse=True)
#        # If alignments exist after filtering report features using alignment pair with largest 
#        # internal segment i.e. first element in sorted list.
#        if alignments:
#            if verbose:
#                [print(x) for x in alignments]
#            if report == 'all':
#                # yield original element
#                yield rec
#            if report in ['split','external']:
#                # yield LTR slice - append "_LTR"
#                extSeg = rec[alignments[0].ref_start:alignments[0].ref_end + 1]
#                extSeg.id = extSeg.id + "_LTR"
#                extSeg.name = extSeg.id
#                extSeg.description = "[" + rec.id + " LTR segment]"
#                yield extSeg
#            if report in ['split','internal']:
#                # yield internal slice - append "_I"
#                intSeg = rec[alignments[0].ref_end:alignments[0].qry_start + 1] 
#                intSeg.id = intSeg.id + "_I"
#                intSeg.name = intSeg.id
#                intSeg.description = "[" + rec.id + " internal segment]"
#                yield intSeg
#        else:
#            # If alignment list empty after filtering print alert and continue
#            print('No LTRs found for candidate element: %s' % rec.id)



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