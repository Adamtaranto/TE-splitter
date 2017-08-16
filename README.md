# TE-exterminate

Extract terminal repeats from retrotransposons (LTRs) or DNA transposons (TIRs). 
Optionally, compose synthetic MITES from complete DNA transposons.  


Currently requires [custom pymummer version](https://github.com/Adamtaranto/pymummer) with wrapper for nucmer option *--diagfactor*.
Pull request pending to update oficial [pymummer](https://github.com/sanger-pathogens/pymummer/pull/28)


# Table of contents

* [Options and usage](#options-and-usage)
    * [Example usage](#example-usage)
    * [Standard options](#standard-options)


# Options and usage  

### Example usage  

For each element in *retroelements.fasta* split into internal and external segments. 
Split segments will be written to *LTR_split_exterminate_output.fasta* with suffix "_I" for internal or "_LTR" for external segments.
LTRs must be at least 10bp in length and share 80% identity and occur within 10bp of each end of the input element

`
./exterminate.py -i retroelements.fasta -p LTR_split --findmode LTR 
`

### Standard options

Run `./exterminate.py --help` to view the program's most commonly used options:

```
usage: ./exterminate.py [-h] -i INFILE [-p PREFIX] [-d OUTDIR]
                        [--findmode {LTR,TIR}]
                        [--splitmode {all,split,internal,external,None}]
                        [--makemites] [-m MAXDIST] [--minid MINID]
                        [--minterm MINTERM] [--minseed MINSEED]
                        [--diagfactor DIAGFACTOR]  


Help:
  -h, --help                        Show this help message and exit


Input:
  -i INFILE, --infile infile        Multifasta containing complete elements. (required)  


Output:
  -p PREFIX, --prefix PREFIX        All output files begin with this string.  (Default:[infile basename])  
  -d OUTDIR, --outdir OUTDIR        Write output files to this directory. (Default: cwd)  


Report settings:
  --findmode                        Type of terminal repeat to identify. (Default: LTR)  
                                      Options: {LTR,TIR}  
  --splitmode                       Options: {all,split,internal,external,None} (Default: split)  
                                      all = Report input sequence as well as internal and external segments.  
                                      split = Report internal and external segments after splitting.  
                                      internal = Report only internal segments.  
                                      external = Report only terminal repeat segments.  
                                      None = Only report synthetic MITES (when --makemites is also set).  
  --makemites                       Attempt to construct synthetic MITE sequences from TIRs.  


Alignment settings:
  --minid MINID                     Minimum identity between terminal repeat pairs. As float. (Default: 80.0)  
  --minterm MINTERM                 Minimum length for a terminal repeat to be considered.  
                                      Equivalent to nucmer "--mincluster" (Default: 10)  
  -m MAXDIST, --maxdist MAXDIST     Terminal repeat candidates must be no more than this many bases from ends of input element. 
                                      Note: Increase this value if you suspect that your element is nested within some flanking sequence. (Default: 10)
  --minseed MINSEED                 Minimum length of a maximal exact match to be included in final match cluster. 
                                      Equivalent to nucmer "--minmatch". (Default: 5)
  --diagfactor DIAGFACTOR           Maximum diagonal difference factor for clustering of matches within nucmer, i.e. diagonal difference / match separation (default 0.20) 
                                      Note: Increase value for greater tolerance of indels between terminal repeats.
```

# License

Software provided under MIT license.


