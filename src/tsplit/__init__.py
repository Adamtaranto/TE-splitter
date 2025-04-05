"""
This is tSplit: A package for analyzing and manipulating transposable elements in genomic sequences.

This package provides tools to identify, extract, and analyze transposable elements (TEs)
in DNA sequences, with a focus on terminal repeats. The main functionality includes:

- Identification of Terminal Inverted Repeats (TIRs) in DNA transposons
- Identification of Long Terminal Repeats (LTRs) in retrotransposons
- Extraction of internal regions and terminal repeats as separate sequences
- Creation of synthetic Miniature Inverted-repeat Transposable Elements (MITEs)
- Support for multiple alignment tools (blastn, nucmer)

The package is built around two main commands:
- tsplit TIR: For processing DNA transposons with Terminal Inverted Repeats
- tsplit LTR: For processing retrotransposons with Long Terminal Repeats

For usage examples, see the documentation at:
https://github.com/adamtaranto/tSplit
"""
