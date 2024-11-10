import os
import logging
import tempfile
import shutil

from Bio import SeqIO
from pymummer import coords_file, nucmer

from tsplit.wrapping import run_cmd, makeBlast
from tsplit.utils import cleanID


def getTIRs(
    fasta_file,
    flankdist=10,
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
        fasta_file (str): Path to the multifasta file containing sequence records.
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

    # Create a unique temporary directory
    tempDir = tempfile.mkdtemp(prefix="tsplit_temp_", dir=temp)
    logging.info(f"Init temp directory: {tempDir}")

    seen_ids = set()

    try:
        # Iterate over each record in the fasta file
        for rec in SeqIO.parse(fasta_file, "fasta"):
            # Log the record name and length
            logging.info(
                f"Processing record {len(seen_ids) + 1}: Name: {rec.id}, Length: {len(rec)}bp"
            )

            # Check for duplicate IDs
            if rec.id in seen_ids:
                logging.error(f"Duplicate record ID found: {rec.id}")
                raise ValueError(f"Duplicate record ID found: {rec.id}")
            seen_ids.add(rec.id)

            # Create temp paths for single element fasta and alignment coords
            tempFasta = os.path.join(tempDir, cleanID(rec.id) + ".fasta")
            tempCoords = os.path.join(
                tempDir, cleanID(rec.id) + "_" + alignTool + ".coords"
            )

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
                run_cmd(cmd, verbose=verbose, workingDir=tempDir)

            # Import coords file to iterator object
            file_reader = coords_file.reader(tempCoords)

            # Exclude hits to self. Also converts iterator output to stable list
            alignments = [hit for hit in file_reader if not hit.is_self_hit()]

            logging.debug(f"NON SELF ALIGNMENTS: {len(alignments)}")

            # Filter hits less than min length (Done internally for nucmer, not blastn.)
            alignments = [
                hit for hit in alignments if hit.ref_end - hit.ref_start >= minterm
            ]

            logging.debug(f"ALIGNMENTS >= minlen {minterm}: {len(alignments)}")

            # Filter for hits on same strand i.e. tandem repeats / LTRs
            alignments = [hit for hit in alignments if not hit.on_same_strand()]

            logging.debug(f"ALIGNMENTS ON OPPOSITE STRANDS: {len(alignments)}")

            # Filter for 5' repeats which begin within x bases of element start
            alignments = [hit for hit in alignments if hit.ref_start <= flankdist]

            logging.debug(
                f"ALIGNMENTS within {flankdist}bp of element start: {len(alignments)}"
            )

            # Scrub overlapping ref / query segments, and also complementary
            # 3' to 5' flank hits
            alignments = [hit for hit in alignments if hit.ref_end < hit.qry_end]
            logging.debug(f"NON-OVERLAPPING ALIGNMENTS: {len(alignments)}")

            # Sort largest to smallest dist between end of ref (subject) and start
            # of query (hit)
            # x.qry_end - x.ref_end =
            # 5'end of right TIR - 3' end of left TIR = length of internal segment
            # TIR pair with largest internal segment (outermost TIRs) is first in list.
            alignments = sorted(
                alignments, key=lambda x: (x.qry_start - x.ref_end), reverse=True
            )
            # If alignments exist after filtering, report features using alignment
            # pair with largest internal segment i.e. first element in sorted list.
            if alignments:
                if verbose:
                    logging.info(f"Alignments found for candidate element: {rec.id}")
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
        if not keeptemp:
            shutil.rmtree(tempDir)
            logging.info(f"Temporary directory deleted: {tempDir}")
        else:
            logging.info(f"Temporary directory retained: {tempDir}")


def getLTRs(
    fasta_file,
    flankdist=10,
    minid=80,
    minterm=10,
    minseed=5,
    diagfactor=0.3,
    report="split",
    temp=None,
    keeptemp=False,
    alignTool="nucmer",
    verbose=True,
):
    """
    Align elements to self and attempt to identify LTRs.

    Args:
        fasta_file (str): Path to the multifasta file containing sequence records.
        flankdist (int): Maximum distance from element start for LTR candidates.
        minid (float): Minimum identity between terminal repeat pairs.
        minterm (int): Minimum length for a terminal repeat to be considered.
        minseed (int): Minimum seed length for nucmer.
        diagfactor (float): Diagonal factor for nucmer.
        report (str): Reporting mode for LTRs.
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

    # Create a unique temporary directory
    tempDir = tempfile.mkdtemp(prefix="tsplit_temp_", dir=temp)
    logging.info(f"Temporary directory created: {tempDir}")

    seen_ids = set()

    try:
        # Iterate over each record in the fasta file
        for rec in SeqIO.parse(fasta_file, "fasta"):
            # Log the record name and length
            logging.info(f"Processing record: {rec.id}, Length: {len(rec)}")

            # Check for duplicate IDs
            if rec.id in seen_ids:
                logging.error(f"Duplicate record ID found: {rec.id}")
                raise ValueError(f"Duplicate record ID found: {rec.id}")

            seen_ids.add(rec.id)

            # Create temp paths for single element fasta and alignment coords
            tempFasta = os.path.join(tempDir, cleanID(rec.id) + ".fasta")
            tempCoords = os.path.join(
                tempDir, cleanID(rec.id) + "_" + alignTool + ".coords"
            )

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
                run_cmd(cmd, verbose=verbose, workingDir=tempDir)

            # Import coords file to iterator object
            file_reader = coords_file.reader(tempCoords)

            # Exclude hits to self. Also converts iterator output to stable list
            alignments = [hit for hit in file_reader if not hit.is_self_hit()]

            # Filter hits less than min length (Done internally for nucmer, not blastn.)
            alignments = [
                hit for hit in alignments if hit.ref_end - hit.ref_start >= minterm
            ]

            # Filter for hits on same strand i.e. tandem repeats / LTRs
            alignments = [hit for hit in alignments if hit.on_same_strand()]

            # Filter for 5' repeats which begin within x bases of element start
            alignments = [hit for hit in alignments if hit.ref_start <= flankdist]

            # Filter for 5' repeats whose 3' match ends within x bases of element end
            alignments = [
                hit for hit in alignments if len(rec) - hit.qry_end <= flankdist
            ]

            # Scrub overlappying ref / query segments, and also complementary 3' to 5' flank hits
            alignments = [hit for hit in alignments if hit.ref_end < hit.qry_start]

            # Sort largest to smallest dist between end of ref (subject) and start of query (hit)
            # x.qry_start (3') - x.ref_end (5') = Length of internal segment
            alignments = sorted(
                alignments, key=lambda x: (x.qry_start - x.ref_end), reverse=True
            )

            # If alignments exist after filtering report features using alignment pair with largest
            # internal segment i.e. first element in sorted list.
            if alignments:
                if verbose:
                    logging.info(f"Alignments found for candidate element: {rec.id}")
                    [print(x) for x in alignments]
                if report == "all":
                    # yield original element
                    yield rec
                if report in ["split", "external"]:
                    # yield LTR slice - append "_LTR"
                    extSeg = rec[alignments[0].ref_start : alignments[0].ref_end + 1]
                    extSeg.id = extSeg.id + "_LTR"
                    extSeg.name = extSeg.id
                    extSeg.description = "[" + rec.id + " LTR segment]"
                    yield extSeg
                if report in ["split", "internal"]:
                    # yield internal slice - append "_I"
                    intSeg = rec[alignments[0].ref_end : alignments[0].qry_start + 1]
                    intSeg.id = intSeg.id + "_I"
                    intSeg.name = intSeg.id
                    intSeg.description = "[" + rec.id + " internal segment]"
                    yield intSeg
            else:
                # If alignment list empty after filtering print alert and continue
                logging.info(f"No LTRs found for candidate element: {rec.id}")
    finally:
        # Clean up the temporary directory if keeptemp is False
        if not keeptemp:
            shutil.rmtree(tempDir)
            logging.info(f"Temporary directory deleted: {tempDir}")
        else:
            logging.info(f"Temporary directory retained: {tempDir}")


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
