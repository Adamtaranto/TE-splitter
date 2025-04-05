# CLI Tutorial

## Overview

The `tsplit` command is used to split transposable elements into their internal and external segments. It can be used for both long terminal repeats (LTRs) and terminal inverted repeats (TIRs).

The command line tools use sequence alignment to identify the boundaries of the segments based on user-defined parameters.

## Example usage

tSplit can be run in two modes: `tsplit LTR` and `tsplit TIR`, for extracting long terminal repeats or terminal inverted repeats, respectively.

Options are the same for each.

## tsplit TIR

For each element in *TIR_element.fa* split into internal and external (TIR) segments.

Split segments will be written to *TIR_split_tsplit_output.fasta* with suffix "_I" for internal or "_TIR" for external segments.

TIRs must be at least 10bp in length and share 80%
identity and occur within 10bp of each end of the input element.

Additionally, synthetic MITEs will be constructed by concatenation of left and right TIRs, with internal segments excised.

```bash
tsplit TIR -i tests/data/TIR_element.fa -p TIR_split --makemites --keeptemp

# Equivalet to defaults
tsplit TIR -i tests/data/TIR_element.fa -p TIR_split --maxdist 10 --minid 80.0 --minterm 10 --method blastn --splitmode split --makemites --keeptemp
```

Output: `TIR_split_tsplit_output.fasta`

## tsplit LTR

For each element in *LTR_retrotransposon.fa* split into internal and external segments.

Split segments will be written to *LTR_split_tsplit_output.fasta* with suffix "_I" for internal or "_LTR" for external segments.

LTRs must be at least 10bp in length and share 80% identity and occur within 10bp of each end of the input element.

```bash
tsplit LTR -i tests/data/LTR_retrotransposon.fa -p LTR_split
```

Output: LTR_split_tsplit_output.fasta
