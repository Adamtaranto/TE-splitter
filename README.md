# tSplit the TE-splitter

Extract terminal repeats from retrotransposons (LTRs) or DNA transposons (TIRs). 
Optionally, compose synthetic MITES from complete DNA transposons.  

# Table of contents

* [Algorithm overview](#algorithm-overview)
* [Options and usage](#options-and-usage)
    * [Installing tSplit](#installing-tsplit)
    * [Example usage](#example-usage)
    * [Standard options](#standard-options)

# Algorithm overview

tSplit attempts to identify terminal repeats in transposable elements by 
first aligning each element to itself using nucmer, and then applying a set of 
tuneable heuristics to select an alignment pair most likely to represent an LTR or TIR.  

  1. Exclude all diagonal/self-matches 
  2. If tsplit-LTR: Retain only alignment pairs on the same strand (tandem repeats)
  3. If tsplit-TIR: Retain only alignment pairs on opposite strands (inverse repeats)
  4. Retain pairs for which the 5' match begins within x bases of element start
     and whose 3' match ends within x bases of element end
  5. Exclude alignment pairs which overlap (potential SSRs)
  6. If multiple candidates remain select alignment pair with largest internal segment 
  (i.e. closest to element ends)

# Options and usage  

### Installing tSplit

Requirements: 
  * [pymummer](https://pypi.python.org/pypi/pymummer) version >= 0.10.3 with wrapper for nucmer option *--diagfactor*.
  * [MUMmer](http://mummer.sourceforge.net/)
  * [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Optional)

Installation options:

```bash
# Install from PyPi:
pip install tsplit

# Clone and install from this repository:
git clone https://github.com/Adamtaranto/TE-splitter.git && cd TE-splitter && pip install -e .
```

### Example usage  

tSplit contains two programs: tsplit-LTR and tsplit-TIR, for extracting long terminal 
repeats and terminal inverted repeats, respectively. Options are the same
for each.  

### tsplit-LTR 

For each element in *retroelements.fasta* split into internal and external segments. 
Split segments will be written to *LTR_split_TE-splitter_output.fasta* with suffix "_I" 
for internal or "_LTR" for external segments. LTRs must be at least 10bp in length and 
share 80% identity and occur within 10bp of each end of the input element.

```bash
tsplit-LTR -i retroelements.fasta -p LTR_split
```

### tsplit-TIR

For each element in *dna-transposons.fasta* split into internal and external (TIR) segments. 
Split segments will be written to *TIR_split_TE-splitter_output.fasta* with suffix "_I" for 
internal or "_TIR" for external segments. TIRs must be at least 10bp in length and share 80% 
identity and occur within 10bp of each end of the input element. Additionally, synthetic 
MITEs will be constructed by concatenation of left and right TIRs, with internal segments 
excised.

```bash
tsplit-TIR -i dna-transposons.fasta -p TIR_split --makemites
```

### Standard options

Run `tsplit-LTR --help` or `tsplit-TIR --help` to view the programs' most commonly used 
options:

```
Usage: tsplit-[LTR or TIR] [-h] -i INFILE [-p PREFIX] [-d OUTDIR]
                        [--splitmode {all,split,internal,external,None}]
                        [--makemites] [--keeptemp] [-v] [-m MAXDIST]
                        [--minid MINID] [--minterm MINTERM] [--minseed MINSEED]
                        [--diagfactor DIAGFACTOR] [--method {blastn,nucmer}]

Help:
  -h, --help         Show this help message and exit.

Input:
  -i, --infile       Multifasta containing complete elements. 
                       (Required)  

Output:
  -p, --prefix       All output files begin with this string.  (Default:[infile basename])  
  -d, --outdir       Write output files to this directory. (Default: cwd)  
  --keeptemp         If set do not remove temp directory on completion.
  -v, --verbose      If set, report progress.

Report settings:
  --splitmode        Options: {all,split,internal,external,None} 
                       all = Report input sequence as well as internal and external segments.  
                       split = Report internal and external segments after splitting.  
                       internal = Report only internal segments.  
                       external = Report only terminal repeat segments.  
                       None = Only report synthetic MITES (when --makemites is also set).  
                       (Default: split)  
  --makemites        Attempt to construct synthetic MITE sequences from TIRs by concatenating 
                       5' and 3' TIRs. Available only in 'tsplit-TIR' mode 

Alignment settings:
  --method          Select alignment tool. Note: blastn may perform better on very short high-identity TRs,
                      while nucmer is more robust to small indels.
                      Options: {blastn,nucmer} 
                      (Default: nucmer)
  --minid           Minimum identity between terminal repeat pairs. As float. 
                      (Default: 80.0)  
  --minterm         Minimum length for a terminal repeat to be considered.  
                      Equivalent to nucmer "--mincluster" 
                      (Default: 10)  
  -m, --maxdist     Terminal repeat candidates must be no more than this many bases from ends of an input element. 
                      Note: Increase this value if you suspect that your element is nested within some flanking sequence. 
                      (Default: 10)
  --minseed         Minimum length of a maximal exact match to be included in final match cluster. 
                      Equivalent to nucmer "--minmatch". 
                      (Default: 5)
  --diagfactor      Maximum diagonal difference factor for clustering of matches within nucmer, 
                      i.e. diagonal difference / match separation 
                      (default 0.20) 
                      Note: Increase value for greater tolerance of indels between terminal repeats.
```

# License

Software provided under MIT license.


