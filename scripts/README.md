**tsd-tir-interval-finder.py**

Takes a fasta file containing one sequence and range coords for two non-overlapping windows in the sequence.

Extracts all kmers from the left window and then searches for matches in the reverse complement of the right window within max 'n' mismatches.

Screen all pairs for direct repeats (posible Target Site Duplications) flanking the left and right TIRs. 

*data/mintest.fa*
```
>tsd_5_TCAGT_k8_AGTCCCGG_int_1-21_33-50
GGTCAGTAGTCCCGGATGATACTATCCATATGCCGGGACTTCAGTGCGC
```

Example:
```bash
python tsd-tir-interval-finder.py --infile data/mintest.fa --left 1-21 --right 33-50 \
-k 8 -n 1 --max-tsd 6 --min-tsd 2 --tsd-mismatches 1 --output test_results.tsv
```

Output:
```
left_pos	left_kmer	right_pos	right_kmer	mismatch_count	left_tsd	right_tsd	tsd_mismatch_count	tsd_len
8	AGTCCCGG	40	AGTCCCGG	0	TCAGT	TCAGT	0	5
```

Output coords are 1-based. Left_pos is begining of left TIR, right_pos is at the outer end of the right TIR (i.e. begining on reverse strand).
