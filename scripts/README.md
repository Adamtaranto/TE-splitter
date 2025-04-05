**tsd-tir-interval-finder.py**

This script is designed to screen two genomic windows for TIR candidates with TSDs using kmer matching.

It takes as input a fasta file containing one sequence and range coords for two non-overlapping windows in the sequence.

Algo overview:

- Extract all kmers from the left window and then searches for matches in the reverse complement of the right window within max 'n' mismatches.

- Screen all pairs for direct repeats (posible Target Site Duplications) flanking the left and right TIRs.

Caveats:

- This tool is not optimised for speed or large search regions.
- Does not perform alignment of candidate TIRs (i.e. only tolerates mismatches, not gaps).
- Only detects terminal k bases of a TIR, actual TIRs may be longer on manual inspection.

*data/mintest.fa*
```
>test_seq   tsd_5_TCAGT_k8_AGTCCCGG_int_1-21_33-50
GGTCAGTAGTCCCGGATGATACTATCCATATGCCGGGACTTCAGTGCGCA
```

Example:
```bash
python tsd-tir-interval-finder.py --infile data/mintest.fa --left 1-21 --right 33-50 \
-k 8 -n 0 --max-tsd 6 --min-tsd 2 --tsd-mismatches 1 --output test_results.tsv --gff test_annotations.gff3
```

Output:
*test_results.tsv*
```
left_pos	left_kmer	right_pos	right_kmer	mismatch_count	left_tsd	right_tsd	tsd_mismatch_count	tsd_len
8	AGTCCCGG	40	AGTCCCGG	0	TCAGT	TCAGT	0	5

```

*test_annotations.gff3*
```
##gff-version 3
test_seq	TIR_finder	search_region	1	21	.	.	.	ID=search_region_L;Name=Left_search_interval
test_seq	TIR_finder	search_region	33	50	.	.	.	ID=search_region_R;Name=Right_search_interval
test_seq	TIR_finder	TIR	8	15	.	+	.	ID=TIR_pair_1_TIR_L;Parent=TIR_pair_1;mismatch_count=0;seq=AGTCCCGG;len=8
test_seq	TIR_finder	TIR	33	40	.	-	.	ID=TIR_pair_1_TIR_R;Parent=TIR_pair_1;mismatch_count=0;seq=AGTCCCGG;len=8
test_seq	TIR_finder	TSD	3	7	.	+	.	ID=TIR_pair_1_TSD_L;Parent=TIR_pair_1_TIR_L;mismatch_count=0;seq=TCAGT;len=5
test_seq	TIR_finder	TSD	41	45	.	+	.	ID=TIR_pair_1_TSD_R;Parent=TIR_pair_1_TIR_R;mismatch_count=0;seq=TCAGT;len=5


```

Output coords are 1-based. Left_pos is begining of left TIR, right_pos is at the outer end of the right TIR (i.e. begining on reverse strand).
