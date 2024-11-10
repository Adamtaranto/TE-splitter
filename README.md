# tSplit the TE-splitter

Extract terminal repeats from retrotransposons (LTRs) or DNA transposons (TIRs). Returns compontent segments of the element for use with transposon mapping tools.

Optionally, `tsplit TIR` can also compose synthetic MITES from complete DNA transposons.  

# Table of contents

* [Algorithm overview](#algorithm-overview)
* [Options and usage](#options-and-usage)
    * [Installing tSplit](#installing-tsplit)
    * [Example usage](#example-usage)
    * [Standard options](#standard-options)

# Algorithm overview

tSplit attempts to identify terminal repeats in transposable elements by 
first aligning each element to itself using `blastn` or `nucmer`, and then applying a set of 
tuneable heuristics to select an alignment pair most likely to represent an LTR or TIR, as follows:

  1. Exclude all diagonal/self-matches 
  2. If `tsplit LTR`: Retain only alignment pairs on the same strand (tandem repeats)
  3. If `tsplit TIR`: Retain only alignment pairs on opposite strands (inverse repeats)
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

tSplit can be run in two modes: `tsplit LTR` and `tsplit TIR`, for extracting long terminal repeats or terminal inverted repeats, respectively. 

Options are the same for each.  

### tsplit-LTR 

For each element in *retroelements.fasta* split into internal and external segments. 

Split segments will be written to *LTR_split_tsplit_output.fasta* with suffix "_I" 
for internal or "_LTR" for external segments. LTRs must be at least 10bp in length and 
share 80% identity and occur within 10bp of each end of the input element.

```bash
tsplit LTR -i retroelements.fasta -p LTR_split
```

### tsplit TIR

For each element in *sample_TIR_element.fa* split into internal and external (TIR) segments. 

Split segments will be written to *TIR_split_tsplit_output.fasta* with suffix "_I" for 
internal or "_TIR" for external segments. TIRs must be at least 10bp in length and share 80% 
identity and occur within 10bp of each end of the input element. 

Additionally, synthetic 
MITEs will be constructed by concatenation of left and right TIRs, with internal segments 
excised.

```bash
tsplit TIR -i tests/data/sample_TIR_element.fa -p TIR_split --makemites --keeptemp
```

# License

Software provided under MIT license.


