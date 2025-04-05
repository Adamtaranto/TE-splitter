#!/usr/bin/env python3

import argparse
import sys
from typing import Dict, Generator, List, NamedTuple, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm


def parse_args() -> argparse.Namespace:
    """Parse and validate command line arguments"""
    parser = argparse.ArgumentParser(
        description='Find inverted repeats and TSDs in DNA sequences'
    )
    parser.add_argument(
        '--infile', required=True, help='Path to input FASTA file containing sequence'
    )
    parser.add_argument(
        '--left', required=True, help='Left genomic interval (format: start-end)'
    )
    parser.add_argument(
        '--right', required=True, help='Right genomic interval (format: start-end)'
    )
    parser.add_argument(
        '-k', type=int, required=True, help='Kmer length for inverted repeat search'
    )
    parser.add_argument(
        '-n',
        type=int,
        default=0,
        help='Number of mismatches to tolerate for inverted repeats (default: 0)',
    )
    parser.add_argument(
        '--max-tsd',
        type=int,
        default=20,
        help='Maximum length of TSD to look for (default: 20)',
    )
    parser.add_argument(
        '--min-tsd',
        type=int,
        default=2,
        help='Minimum length of TSD to report (default: 2)',
    )
    parser.add_argument(
        '--tsd-mismatches',
        type=int,
        default=1,
        help='Number of mismatches to tolerate in TSDs (default: 1)',
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Output file path for tab-delimited results (optional)',
    )
    parser.add_argument(
        '--gff', type=str, help='Output file path for GFF3 format results (optional)'
    )

    return parser.parse_args()


class GenomicInterval(NamedTuple):
    """Store start and end coordinates for a genomic interval"""

    start: int  # 1-based start position
    end: int  # 1-based end position


def load_sequence(fasta_path: str, min_length: int) -> Seq:
    """
    Load a single sequence from a FASTA file and validate it.

    Args:
        fasta_path: Path to FASTA file
        min_length: Minimum required sequence length

    Returns:
        Bio.Seq.Seq: The loaded sequence

    Raises:
        ValueError: If file contains multiple sequences or sequence is too short
        FileNotFoundError: If file cannot be opened
    """
    try:
        sequences = list(SeqIO.parse(fasta_path, 'fasta'))

        if not sequences:
            raise ValueError(f'No sequences found in {fasta_path}')

        if len(sequences) > 1:
            raise ValueError(
                f'Multiple sequences found in {fasta_path}. Expected single sequence.'
            )

        seq = sequences[0].seq

        if len(seq) < min_length:
            raise ValueError(
                f'Sequence in {fasta_path} is too short. '
                f'Length: {len(seq)}, Required: {min_length}'
            )

        return seq

    except FileNotFoundError:
        raise FileNotFoundError(f'Could not open FASTA file: {fasta_path}') from None


def parse_interval(
    interval_str: str, seq_length: int, interval_name: str
) -> GenomicInterval:
    """
    Parse a string of format 'start-end' into genomic coordinates.

    Args:
        interval_str: String in format 'start-end'
        seq_length: Length of the reference sequence
        interval_name: Name of interval for error messages

    Returns:
        GenomicInterval with start and end coordinates

    Raises:
        ValueError: If format is invalid or coordinates are out of bounds
    """
    try:
        start, end = map(int, interval_str.split('-'))
    except ValueError:
        raise ValueError(
            f'Invalid {interval_name} interval format: {interval_str}. '
            'Expected format: start-end'
        ) from None

    if start > end:
        raise ValueError(
            f'Invalid {interval_name} interval: start ({start}) > end ({end})'
        )

    if end > seq_length:
        raise ValueError(
            f'{interval_name.capitalize()} interval end ({end}) exceeds '
            f'sequence length ({seq_length})'
        )

    if start < 1:
        raise ValueError(
            f'{interval_name.capitalize()} interval start ({start}) must be ≥ 1'
        )

    return GenomicInterval(start, end)


def check_interval_overlap(left: GenomicInterval, right: GenomicInterval) -> bool:
    """Check if two genomic intervals overlap"""
    return (left.start <= right.end) and (right.start <= left.end)


def extract_sequence_blocks(
    sequence: Seq, left: GenomicInterval, right: GenomicInterval
) -> Tuple[Seq, Seq]:
    """
    Extract left and right sequence blocks from genomic intervals.
    Convert from 1-based to 0-based coordinates for slicing.

    Args:
        sequence: Full DNA sequence
        left: Left genomic interval (1-based coords)
        right: Right genomic interval (1-based coords)

    Returns:
        Tuple of (left_block, right_block)
    """
    # Convert to 0-based coordinates for slicing
    print('\nExtracting sequence blocks...')
    print(f'Left interval: {left.start}-{left.end}')
    print(f'Right interval: {right.start}-{right.end}')

    left_block = sequence[left.start - 1 : left.end]
    right_block = sequence[right.start - 1 : right.end]

    return left_block, right_block


def calculate_hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculate the Hamming distance between two sequences.
    'N' in either sequence always counts as a mismatch.

    Args:
        seq1: First sequence
        seq2: Second sequence

    Returns:
        int: Number of positions at which sequences differ

    Raises:
        ValueError: If sequences have different lengths
    """
    if len(seq1) != len(seq2):
        raise ValueError('Sequences must have equal length')

    return sum(1 for a, b in zip(seq1, seq2) if a != b or 'N' in (a, b))


def find_inverted_repeats(
    left_block: Seq, right_block: Seq, k: int, n_mismatches: int
) -> Generator[Tuple[str, int, str, int], None, None]:
    """
    Find inverted repeats between two DNA sequence blocks allowing for mismatches.

    Args:
        left_block: Left DNA sequence block (Seq object)
        right_block: Right DNA sequence block (Seq object)
        k: Length of kmers to use
        n_mismatches: Number of mismatches to tolerate

    Yields:
        Tuples of (left_kmer, left_position, right_kmer, right_position)
    """
    # Reverse complement right block
    right_seq = right_block.reverse_complement()

    # Build left kmer index
    left_index: Dict[str, List[int]] = {}
    print('Indexing left block kmers...')
    for i in tqdm(range(len(left_block) - k + 1)):
        kmer = str(left_block[i : i + k])
        if 'N' not in kmer:
            if kmer not in left_index:
                left_index[kmer] = []
            left_index[kmer].append(i)

    # Build right kmer index
    right_index: Dict[str, List[int]] = {}
    print('Indexing right block kmers...')
    for i in tqdm(range(len(right_seq) - k + 1)):
        kmer = str(right_seq[i : i + k])
        if 'N' not in kmer:
            if kmer not in right_index:
                right_index[kmer] = []
            right_index[kmer].append(i)

    # Find mismatch pairs
    mismatch_pairs: Dict[str, List[str]] = {}
    print(f'Finding inverted repeat pairs within {n_mismatches} mismatches...')
    for left_kmer in tqdm(left_index.keys()):
        mismatch_pairs[left_kmer] = []
        for right_kmer in right_index.keys():
            if calculate_hamming_distance(left_kmer, right_kmer) <= n_mismatches:
                mismatch_pairs[left_kmer].append(right_kmer)

    # Generate all valid kmer position pairs
    print('Generating pair index positions...')
    for left_kmer in tqdm(mismatch_pairs.keys()):
        left_positions = left_index[left_kmer]
        for right_kmer in mismatch_pairs[left_kmer]:
            right_positions = right_index[right_kmer]
            for left_pos in left_positions:
                for right_pos in right_positions:
                    yield (left_kmer, left_pos, right_kmer, right_pos)


def find_tsds(
    repeat_pairs: List[Tuple[str, int, str, int]],
    left_block: Seq,
    right_block: Seq,
    max_tsd_len: int,
    n_mismatches: int,
    min_tsd_len: int,
) -> Generator[Tuple[str, int, str, int, str, str], None, None]:
    """
    Find Target Site Duplications (TSDs) in the flanking sequences of inverted repeats.

    Args:
        repeat_pairs: List of tuples from find_inverted_repeats()
        left_block: Left DNA sequence block (Seq object)
        right_block: Right DNA sequence block (Seq object)
        max_tsd_len: Maximum length of TSD to look for
        n_mismatches: Number of mismatches to tolerate
        min_tsd_len: Minimum length of TSD to report

    Yields:
        Tuples of (left_kmer, left_position, right_kmer, right_position, left_tsd, right_tsd)
        Note: right_tsd is reverse complemented in the output
    """
    # Convert to Biopython Seq objects and reverse complement right block
    right_seq = right_block.reverse_complement()

    print('\nSearching for TSDs...')
    for repeat in tqdm(repeat_pairs):
        left_kmer, left_pos, right_kmer, right_pos = repeat

        # Try decreasing lengths of flanking sequence
        for k in range(max_tsd_len, min_tsd_len - 1, -1):
            # Check if we have enough sequence to extract flanks
            if left_pos < k or right_pos < k:
                continue

            # Extract flanking sequences
            left_flank = str(left_block[left_pos - k : left_pos])
            right_flank = str(right_seq[right_pos - k : right_pos])

            # Reverse complement the right flank for comparison
            right_flank_rc = str(Seq(right_flank).reverse_complement())

            # Calculate hamming distance
            try:
                distance = calculate_hamming_distance(left_flank, right_flank_rc)

                # If we found a match within tolerance
                if distance <= n_mismatches:
                    yield (
                        left_kmer,
                        left_pos,
                        right_kmer,
                        right_pos,
                        left_flank,
                        right_flank_rc,
                    )  # Return reverse complemented right_flank
                    break  # Stop looking for shorter TSDs once we find a match

            except ValueError:
                continue  # Skip if sequences have different lengths


def correct_reverse_complement_coordinate(position: int, sequence_length: int) -> int:
    """
    Convert a coordinate from a reverse complemented sequence back to its original position.

    Args:
        position: Position in the reverse complemented sequence (0-based)
        sequence_length: Length of the original sequence

    Returns:
        int: Corrected position in the original sequence
    """
    return sequence_length - position


def correct_output_coordinates(
    results: List[Tuple[str, int, str, int, str, str]],
    left_start: int,
    right_start: int,
    right_seq_length: int,
) -> List[Tuple[str, int, str, int, str, str]]:
    """
    Correct output coordinates by adding the original interval start positions.

    Args:
        results: List of result tuples from find_tsds()
        left_start: Start coordinate of left interval (1-based)
        right_start: Start coordinate of right interval (1-based)
        right_seq_length: Length of right sequence block

    Returns:
        List of corrected result tuples
    """
    corrected_results = []

    for result in results:
        left_kmer, left_pos, right_kmer, right_pos, left_tsd, right_tsd = result

        # Correct left position: add left interval start (convert to 1-based)
        corrected_left_pos = left_pos + left_start

        # Correct right position: add right interval start (convert to 1-based)
        corrected_right_pos = correct_reverse_complement_coordinate(
            right_pos, right_seq_length
        ) + (right_start - 1)

        corrected_results.append(
            (
                left_kmer,
                corrected_left_pos,
                right_kmer,
                corrected_right_pos,
                left_tsd,
                right_tsd,
            )
        )

    return corrected_results


def write_results(
    results: List[Tuple[str, int, str, int, str, str]],
    output_path: str,
    right_seq_length: int,
) -> None:
    """
    Write results to a tab-delimited file.

    Args:
        results: List of result tuples from find_tsds()
        output_path: Path to output file
        right_seq_length: Length of right sequence for coordinate correction

    Format:
        left_pos left_kmer right_pos right_kmer mismatch_count left_tsd right_tsd
        tsd_mismatch_count tsd_len

    Results are sorted by:
        1. tsd_len (descending)
        2. left_pos (ascending)
    """
    print('\nProcessing results...')

    # Calculate mismatches and prepare sorted results
    processed_results = []
    for result in tqdm(results, desc='Calculating mismatch counts'):
        left_kmer, left_pos, right_kmer, right_pos, left_tsd, right_tsd = result

        # Calculate mismatch counts
        kmer_mismatches = calculate_hamming_distance(left_kmer, right_kmer)
        tsd_mismatches = calculate_hamming_distance(left_tsd, right_tsd)

        processed_results.append(
            {
                'left_pos': left_pos,
                'left_kmer': left_kmer,
                'right_pos': right_pos,
                'right_kmer': right_kmer,
                'mismatch_count': kmer_mismatches,
                'left_tsd': left_tsd,
                'right_tsd': right_tsd,
                'tsd_mismatch_count': tsd_mismatches,
                'tsd_len': len(left_tsd),
            }
        )

    # Sort results
    print('\nSorting results...')
    sorted_results = sorted(
        processed_results, key=lambda x: (-x['tsd_len'], x['left_pos'])
    )

    # Write to file
    print(f'\nWriting results to {output_path}...')
    with open(output_path, 'w') as f:
        # Write header
        header = [
            'left_pos',
            'left_kmer',
            'right_pos',
            'right_kmer',
            'mismatch_count',
            'left_tsd',
            'right_tsd',
            'tsd_mismatch_count',
            'tsd_len',
        ]
        f.write('\t'.join(header) + '\n')

        # Write sorted results
        for result in tqdm(sorted_results, desc='Writing output'):
            line = [
                str(result['left_pos']),
                result['left_kmer'],
                str(result['right_pos']),
                result['right_kmer'],
                str(result['mismatch_count']),
                result['left_tsd'],
                result['right_tsd'],
                str(result['tsd_mismatch_count']),
                str(result['tsd_len']),
            ]
            f.write('\t'.join(line) + '\n')

    print('Done writing results.')


def write_gff(
    results: List[dict],
    output_path: str,
    k: int,
    sequence_id: str,
    left_interval: str,
    right_interval: str,
) -> None:
    """
    Write TIR and TSD pairs to GFF3 format, including search interval annotations.

    Args:
        results: List of processed result dictionaries
        output_path: Path to output GFF file
        k: Length of TIRs (kmer length)
        sequence_id: ID of the reference sequence
        left_interval: Left interval string (format: start-end)
        right_interval: Right interval string (format: start-end)

    GFF3 Format:
        seqid, source, type, start, end, score, strand, phase, attributes
    """
    print(f'\nWriting GFF to {output_path}...')

    # Parse interval coordinates
    left_start, left_end = map(int, left_interval.split('-'))
    right_start, right_end = map(int, right_interval.split('-'))

    # Calculate total number of records for progress bar
    total_records = (
        len(results) * 4 + 2
    )  # 4 records per result (2 TIRs + 2 TSDs) + 2 intervals

    with open(output_path, 'w') as f:
        # Write GFF header
        f.write('##gff-version 3\n')

        # Initialize progress bar for all records
        with tqdm(total=total_records, desc='Writing GFF records') as pbar:
            # Write search interval annotations
            f.write(
                f'{sequence_id}\tTIR_finder\tsearch_region\t{left_start}\t{left_end}\t'
                f'.\t.\t.\tID=search_region_L;Name=Left_search_interval\n'
            )
            pbar.update(1)

            f.write(
                f'{sequence_id}\tTIR_finder\tsearch_region\t{right_start}\t{right_end}\t'
                f'.\t.\t.\tID=search_region_R;Name=Right_search_interval\n'
            )
            pbar.update(1)

            # Process each TIR pair
            for i, result in enumerate(results):
                pair_id = f'TIR_pair_{i + 1}'

                # Calculate coordinates (convert to 1-based for GFF)
                left_tir_start = result['left_pos']
                left_tir_end = left_tir_start + k - 1
                right_tir_start = result['right_pos'] - k + 1
                right_tir_end = result['right_pos']

                # Calculate TSD coordinates
                tsd_len = result['tsd_len']
                left_tsd_start = left_tir_start - tsd_len
                left_tsd_end = left_tir_start - 1
                right_tsd_start = right_tir_end + 1
                right_tsd_end = right_tir_end + tsd_len

                # Write left TIR
                left_tir_attrs = (
                    f'ID={pair_id}_TIR_L;'
                    f'Parent={pair_id};'
                    f'mismatch_count={result["mismatch_count"]};'
                    f'seq={result["left_kmer"]};'
                    f'len={k}'
                )
                f.write(
                    f'{sequence_id}\tTIR_finder\tTIR\t{left_tir_start}\t{left_tir_end}\t'
                    f'.\t+\t.\t{left_tir_attrs}\n'
                )
                pbar.update(1)

                # Write right TIR
                right_tir_attrs = (
                    f'ID={pair_id}_TIR_R;'
                    f'Parent={pair_id};'
                    f'mismatch_count={result["mismatch_count"]};'
                    f'seq={result["right_kmer"]};'
                    f'len={k}'
                )
                f.write(
                    f'{sequence_id}\tTIR_finder\tTIR\t{right_tir_start}\t{right_tir_end}\t'
                    f'.\t-\t.\t{right_tir_attrs}\n'
                )
                pbar.update(1)

                # Write left TSD if exists
                if tsd_len > 0:
                    left_tsd_attrs = (
                        f'ID={pair_id}_TSD_L;'
                        f'Parent={pair_id}_TIR_L;'
                        f'mismatch_count={result["tsd_mismatch_count"]};'
                        f'seq={result["left_tsd"]};'
                        f'len={tsd_len}'
                    )
                    f.write(
                        f'{sequence_id}\tTIR_finder\tTSD\t{left_tsd_start}\t{left_tsd_end}\t'
                        f'.\t+\t.\t{left_tsd_attrs}\n'
                    )
                    pbar.update(1)

                    # Write right TSD
                    right_tsd_attrs = (
                        f'ID={pair_id}_TSD_R;'
                        f'Parent={pair_id}_TIR_R;'
                        f'mismatch_count={result["tsd_mismatch_count"]};'
                        f'seq={result["right_tsd"]};'
                        f'len={tsd_len}'
                    )
                    f.write(
                        f'{sequence_id}\tTIR_finder\tTSD\t{right_tsd_start}\t{right_tsd_end}\t'
                        f'.\t+\t.\t{right_tsd_attrs}\n'
                    )
                    pbar.update(1)

    print('Done writing GFF.')


def main():
    """Main function to run the inverted repeat and TSD search"""
    # Parse command line arguments
    args = parse_args()

    # Program called with cmd line arguments
    print('\nProgram called with cmd line arguments:')
    print('cmd: {0}'.format(' '.join(sys.argv)))

    try:
        # Load sequence from FASTA file
        sequences = list(SeqIO.parse(args.infile, 'fasta'))

        if not sequences:
            raise ValueError(f'No sequences found in {args.infile}')

        if len(sequences) > 1:
            raise ValueError(
                f'Multiple sequences found in {args.infile}. Expected single sequence.'
            )

        sequence = sequences[0].seq

        # Parse and validate intervals
        left_interval = parse_interval(args.left, len(sequence), 'left')
        right_interval = parse_interval(args.right, len(sequence), 'right')

        # Check for interval overlap
        if check_interval_overlap(left_interval, right_interval):
            print(
                'Warning: Left and right intervals overlap. This may affect results.',
                file=sys.stderr,
            )

        # Extract sequence blocks
        left_seq, right_seq = extract_sequence_blocks(
            sequence, left_interval, right_interval
        )

        # Validate sequence lengths
        if len(left_seq) < args.k or len(right_seq) < args.k:
            raise ValueError(f'Sequence blocks must be at least {args.k} bases long')

        # Find inverted repeats
        print('\nSearching for inverted repeats...')
        repeat_pairs = list(find_inverted_repeats(left_seq, right_seq, args.k, args.n))

        # Find TSDs
        results = list(
            find_tsds(
                repeat_pairs,
                left_seq,
                right_seq,
                args.max_tsd,
                args.tsd_mismatches,
                args.min_tsd,
            )
        )

        # Correct coordinates to original sequence positions
        corrected_results = correct_output_coordinates(
            results, left_interval.start, right_interval.start, len(right_seq)
        )

        print(f'\nTotal TIR pairs: {len(corrected_results)}')

        # Write results to file if output path is provided
        if args.output:
            write_results(corrected_results, args.output, len(right_seq))

        # Write GFF if requested
        if args.gff:
            # Get sequence ID from FASTA header
            sequence_id = sequences[0].id

            # Convert results to dictionary format for GFF writing
            processed_results = []
            for result in corrected_results:
                left_kmer, left_pos, right_kmer, right_pos, left_tsd, right_tsd = result
                processed_results.append(
                    {
                        'left_pos': left_pos,
                        'right_pos': right_pos,
                        'left_kmer': left_kmer,
                        'right_kmer': right_kmer,
                        'left_tsd': left_tsd,
                        'right_tsd': right_tsd,
                        'mismatch_count': calculate_hamming_distance(
                            left_kmer, right_kmer
                        ),
                        'tsd_mismatch_count': calculate_hamming_distance(
                            left_tsd, right_tsd
                        ),
                        'tsd_len': len(left_tsd),
                    }
                )

            write_gff(
                processed_results, args.gff, args.k, sequence_id, args.left, args.right
            )

        # Print results to console
        # count = 0
        # for result in corrected_results:
        #    left_kmer, left_pos, right_kmer, right_pos, left_tsd, right_tsd = result

        #    count += 1

        #    print(
        #        f"TIR #{count}:  {left_kmer}/{right_kmer} at positions {left_pos},{right_pos}"
        #    )

        #    print(f"TSDs: {left_tsd}/{right_tsd} Length: {len(left_tsd)}")

    except (ValueError, FileNotFoundError) as e:
        print(f'Error: {str(e)}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
