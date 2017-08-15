#!/usr/bin/env python
#python 2.7.5
#exterminate.py
#Version 1.0 Adam Taranto, August 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

####################################################################################
# Extract terminal repeats from retrotransposons (LTRs) or DNA transposons (TIRs). #
# Compose synthetic MITES from complete DNA transposons.                           #
####################################################################################

import os
import re
import sys
import shutil
from datetime import datetime
import argparse
from collections import Counter
from Bio import SeqIO
from pymummer import coords_file, alignment, nucmer

def dochecks(args):
	"""Make outDir if does not exist else set to current dir."""
	if args.outdir:
		absOutDir = os.path.abspath(args.outdir)
		if not os.path.isdir(args.outdir):
			os.makedirs(absOutDir)
		outDir = absOutDir
	else:
		outDir = os.getcwd() 

	#Make temp directory
	tempDir = os.path.join(os.getcwd(),"temp_" + getTimestring())
	os.makedirs(tempDir)

	return outDir,tempDir

def getTimestring():
	"""..."""
	(dt, micro) = datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.')
	dt = "%s%03d" % (dt, int(micro) / 1000)
	return dt

def checkUniqueID(records):
	"""..."""
	seqIDs = [records[x].id for x in range(len(records))]
	IDcounts = Counter(seqIDs)
	duplicates = [k for k, v in IDcounts.items() if v > 1]
	if duplicates:
		print("Input sequence IDs not unique. Quiting.")
		print(duplicates)
		sys.exit(1)
	else:
		pass

def manageTemp(record=None,tempPath=None,scrub=False):
	"""..."""
	if scrub and tempPath:
		try:
			os.remove(tempPath)
		except OSError:
			pass
	else:
		with open(tempPath, "w") as f:
			SeqIO.write(record, f, "fasta")

def importFasta(file):
	"""..."""
	# Read in elements from multifasta file, convert seqrecord iterator to list.
	records = list(SeqIO.parse(file, "fasta"))
	# Check names are unique
	checkUniqueID(records)
	# If unique, return record list.
	return records

def cleanID(s):
	"""..."""
	s = re.sub(r"[^\w\s]", '', s)
	s = re.sub(r"\s+", '_', s)
	return s

def getLTRs(elements=None,flankdist=10,minterm=10,minid=90,report='split',temp=None):
	"""..."""
	# Set temp directory to cwd if none.
	if not temp:
		temp = os.getcwd()
	# For each candidate LTR element
	for rec in elements:
		# Create temp paths for single element fasta and alignment coords
		tempFasta = os.path.join(temp, cleanID(rec.id) + '.fasta')
		tempCoords = tempFasta + '.coords'
		# Write current element to single fasta
		manageTemp(record=rec,tempPath=tempFasta,scrub=False)
		# Compose Nucmer script for current element vs self
		runner = nucmer.Runner(tempFasta, tempFasta, tempCoords,
						min_id=minid, min_length=minterm, 
						breaklen=100, maxmatch=True, 
						simplify=False)
		# Execute nucmer
		runner.run()
		# Import coords file to iterator object
		file_reader = coords_file.reader(tempCoords)
		# Exclude hits to self. Also converts iterator output to stable list
		alignments = [hit for hit in file_reader if not hit.is_self_hit()]
		# Filter for hits on same strand i.e. tandem repeats / LTRs
		alignments = [hit for hit in alignments if hit.on_same_strand()]
		# Filter for 5' repeats which begin within x bases of element start
		alignments = [hit for hit in alignments if hit.ref_start <= flankdist]
		# Scrub overlappying ref / query segments, and also complementary 3' to 5' flank hits
		alignments = [hit for hit in alignments if hit.ref_end < hit.qry_start]
		# Sort largest to smallest dist between end of ref (subject) and start of query (hit)
		# x.qry_start - x.ref_end = length of internal segment
		alignments = sorted(alignments, key=lambda x: (x.qry_start - x.ref_end), reverse=True)
		# If alignments exist after filtering report features using alignment pair with largest 
		# internal segment i.e. first element in sorted list.
		if alignments:
			if report == 'all':
				yield rec
			if report in ['split','external']:
				yield rec[alignments[0].ref_start:alignments[0].ref_end] # yield LTR slice - append "_LTR"
			if report in ['split','internal']:
				yield rec[alignments[0].ref_end:alignments[0].qry_start] # yield internal slice - append "_I"
		else:
			# If alignment list empty after filtering print alert and continue
			print('No LTRs found for candidate element: %s' % rec.id)
			pass
		# Scrub single fasta and coords file for current element.
		manageTemp(tempPath=tempFasta,scrub=True)
		manageTemp(tempPath=tempCoords,scrub=True)

def getTIRs(elements=None,flankdist=10,minterm=10,minid=90,mites=False,report='split',temp=None):
	"""..."""
	# Set temp directory to cwd if none.
	if not temp:
		temp = os.getcwd()
	# For each candidate LTR element
	for rec in elements:
		# Create temp paths for single element fasta and alignment coords
		tempFasta = os.path.join(temp, cleanID(rec.id) + '.fasta')
		tempCoords = tempFasta + '.coords'
		# Write current element to single fasta
		manageTemp(record=rec,tempPath=tempFasta,scrub=False)
		# Compose Nucmer script for current element vs self
		runner = nucmer.Runner(tempFasta, tempFasta, tempCoords,
						min_id=minid, min_length=minterm, 
						breaklen=100, maxmatch=True, 
						simplify=False)
		# Execute nucmer
		runner.run()
		# Import coords file to iterator object
		file_reader = coords_file.reader(tempCoords)
		# Exclude hits to self. Also converts iterator output to stable list
		alignments = [hit for hit in file_reader if not hit.is_self_hit()]
		# Filter for hits on same strand i.e. tandem repeats / LTRs
		alignments = [hit for hit in alignments if not hit.on_same_strand()]
		# Filter for 5' repeats which begin within x bases of element start
		alignments = [hit for hit in alignments if hit.ref_start <= flankdist]
		# Scrub overlappying ref / query segments, and also complementary 3' to 5' flank hits
		alignments = [hit for hit in alignments if hit.ref_end < hit.qry_start]
		# Sort largest to smallest dist between end of ref (subject) and start of query (hit)
		# x.qry_start - x.ref_end = length of internal segment
		alignments = sorted(alignments, key=lambda x: (x.qry_start - x.ref_end), reverse=True)
		# If alignments exist after filtering report features using alignment pair with largest 
		# internal segment i.e. first element in sorted list.
		if alignments:
			if report == 'all':
				yield rec
			if report in ['split','external']:
				yield rec[alignments[0].ref_start:alignments[0].ref_end] # yield LTR slice - append "_TIR"
			if report in ['split','internal']:
				yield rec[alignments[0].ref_end:alignments[0].qry_start] # yield internal slice - append "_I"
			if mites:
				# Reassemble TIRs into hypothetical MITEs
				yield str(rec[alignments[0].ref_start:alignments[0].ref_end]) + str(rec[alignments[0].qry_start:alignments[0].qry_end])
		else:
			# If alignment list empty after filtering print alert and continue
			print('No LTRs found for candidate element: %s' % rec.id)
			pass
		# Scrub single fasta and coords file for current element.
		manageTemp(tempPath=tempFasta,scrub=True)
		manageTemp(tempPath=tempCoords,scrub=True)


def segWrite(outhandle,segs=None):
	"""..."""
	if segs:
		with open(outhandle, "w") as f:
			for seq in segs:
				SeqIO.write(seq, 'fasta', f)
	else:
		pass

def mainArgs():
	parser = argparse.ArgumentParser(
			description='Extract terminal repeats from retrotransposons (LTRs) or DNA transposons (TIRs). Compose synthetic MITES from complete DNA transposons.',
			prog='TE-exterminate')
	parser.add_argument('-i', '--infile',
							type=str,
							required=True,
							default=None,
							help='Multifasta containing complete elements.'
							)
	parser.add_argument('-p', '--prefix',
							type=str,
							default='exterminate_split_repeats',
							help='All output files begin with this string.'
							)
	parser.add_argument('-d', '--outdir',
							type=str,
							default=None,
							help='Write output files to this directory.'
							)	
	parser.add_argument('--findmode',
							default='LTR',
							choices=['LTR','TIR'],
							help='Type of terminal repeat to identify. LTR or TIR.')
	parser.add_argument('--splitmode',
							default='split',
							choices=['all','split','internal','flanks',None],
							help='All= Report input sequence as well as internal and external segments. \
							split= Report internal and external segments after splitting. \
							internal = Report only internal segments \
							flanks= Report only flanking segments.')
	parser.add_argument('--makemites',
							action='store_true',
							default=False,
							help='Attempt to construct synthetic MITE sequences from TIRs.')
	parser.add_argument('-m', '--maxdist',
							type=int,
							default=10,
							help='Terminal repeat candidates must be no more than this many bases from end of input element.'
							)
	parser.add_argument('--minid',
							type=float,
							default=95.0,
							help='Minimum identity between terminal repeat pairs. As float. Default: 95.0'
							)
	parser.add_argument('--minterm',
							type=int,
							default=10,
							help='Minimum length for a terminal repeat to be considered.'
							)
	args = parser.parse_args()
	return args

def main(args):
	outDir,tempdir = dochecks(args)

	elements = importFasta(args.infile)

	outfile = args.prefix + "_exterm_extracted.fasta"
	outpath = os.path.join(outDir,outfile)

	if args.findmode == 'LTR':
		segments = getLTRs(elements=elements,flankdist=args.maxdist,minterm=args.minterm,minid=args.minid,report=args.splitmode,temp=tempdir)
		segWrite(outpath,segs=segments)
	elif args.findmode == 'TIR':
		segments = getTIRs(elements=elements,flankdist=args.maxdist,minterm=args.minterm,minid=args.minid,mites=args.makemites,report=args.splitmode,temp=tempdir)
		segWrite(outpath,segs=segments)


	# Remove temp directory
	shutil.rmtree(tempdir)

if __name__== '__main__':
	args = mainArgs()
	main(args)
"""
##Filter for within x% of end
[x.ref_start for x in alignments] #Note: Idexed from 1
[x.ref_end for x in alignments] #Note: Idexed from 1
[x.qry_start for x in alignments] #Note: Idexed from 1
[x.qry_end for x in alignments] #Note: Idexed from 1
[x.hit_length_ref for x in alignments]
[x.hit_length_qry for x in alignments]
[x.percent_identity for x in alignments]
[x.ref_length for x in alignments]
[x.qry_length for x in alignments]
[x.frame for x in alignments]
[x.ref_name for x in alignments]
[x.qry_name for x in alignments]

#coord.reverse_query()
#coord.reverse_reference()
#coord.qry_coords()
#coord.ref_coords()

"""
"""
fields = line.rstrip().split('\t')

"""












