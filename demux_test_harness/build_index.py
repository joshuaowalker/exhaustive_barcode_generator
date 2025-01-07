#!/usr/bin/env python3

"""
Generate Index.txt file for dual-indexed fungal amplicon sequencing with nanopore.
Creates all permutations of forward and reverse barcodes with fixed primers.
"""

import argparse
from pathlib import Path
from typing import List, TextIO
import sys

# Constants for primers
FWD_PRIMER = "CTTGGTCATTTAGAGGAAGTAA"  # ITS1-F
REV_PRIMER = "TCCTCCGCTTATTGATATGC"  # ITS4


def read_barcodes_from_file(file_path: Path) -> List[str]:
    """Read barcodes from a file, one per line."""
    with open(file_path) as f:
        return [line.strip().upper() for line in f if line.strip()]


def write_index_file(fwd_barcodes: List[str], rev_barcodes: List[str],
                     output: TextIO, include_header: bool = True) -> None:
    """Generate and write all barcode permutations to the output file."""

    # Write header if requested
    if include_header:
        output.write("SampleID\tFwIndex\tFwPrimer\tRvIndex\tRvPrimer\n")

    # Generate all permutations
    for i, fwd in enumerate(fwd_barcodes, 1):
        for j, rev in enumerate(rev_barcodes, 1):
            # Format: ONTxx.yy where xx is fwd index and yy is rev index
            sample_id = f"ONT{i:02d}.{j:02d}"

            # Write the line with tab separation
            output.write(f"{sample_id}\t{fwd}\t{FWD_PRIMER}\t{rev}\t{REV_PRIMER}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Generate Index.txt file for dual-indexed nanopore sequencing')

    # Input options group
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--fwd-file', type=Path,
                             help='File containing forward barcodes (one per line)')
    input_group.add_argument('--fwd-barcodes', nargs='+',
                             help='Space-separated list of forward barcodes')

    rev_group = parser.add_mutually_exclusive_group(required=True)
    rev_group.add_argument('--rev-file', type=Path,
                           help='File containing reverse barcodes (one per line)')
    rev_group.add_argument('--rev-barcodes', nargs='+',
                           help='Space-separated list of reverse barcodes')

    # Output options
    parser.add_argument('-o', '--output', type=Path,
                        help='Output file (default: stdout)')
    parser.add_argument('--no-header', action='store_true',
                        help='Omit header line from output')

    args = parser.parse_args()

    # Get forward barcodes
    if args.fwd_file:
        fwd_barcodes = read_barcodes_from_file(args.fwd_file)
    else:
        fwd_barcodes = [bc.upper() for bc in args.fwd_barcodes]

    # Get reverse barcodes
    if args.rev_file:
        rev_barcodes = read_barcodes_from_file(args.rev_file)
    else:
        rev_barcodes = [bc.upper() for bc in args.rev_barcodes]

    # Validate barcodes are not empty
    if not fwd_barcodes or not rev_barcodes:
        parser.error("No barcodes provided")

    # Validate all barcodes contain only valid nucleotides
    valid_bases = set('ACGT')
    for bc in fwd_barcodes + rev_barcodes:
        if not set(bc).issubset(valid_bases):
            parser.error(f"Invalid barcode sequence: {bc}")

    # Validate consistent lengths
    fwd_lens = {len(bc) for bc in fwd_barcodes}
    rev_lens = {len(bc) for bc in rev_barcodes}
    if len(fwd_lens) > 1:
        parser.error("Forward barcodes have inconsistent lengths")
    if len(rev_lens) > 1:
        parser.error("Reverse barcodes have inconsistent lengths")

    # Write output
    if args.output:
        with open(args.output, 'w') as f:
            write_index_file(fwd_barcodes, rev_barcodes, f, not args.no_header)
    else:
        write_index_file(fwd_barcodes, rev_barcodes, sys.stdout, not args.no_header)


if __name__ == "__main__":
    main()