#!/usr/bin/env python3

"""
Nanopore Read Simulator for testing barcode demultiplexing

This program generates synthetic nanopore reads based on an Index.txt file,
simulating errors and structural variations commonly seen in nanopore data.
"""

import argparse
import csv
import random
from typing import List, Tuple, Optional
from dataclasses import dataclass
from Bio.Seq import Seq, reverse_complement
import logging
from pathlib import Path

# Golden string delimiters
START_DELIMITER = "AAAAAAAAAAGGGGGGGGGG"
END_DELIMITER =   "GGGGGGGGGGAAAAAAAAAA"

ERROR_RANGE = 100

@dataclass
class SpecimenInfo:
    """Container for specimen barcode and primer information"""
    id: str
    fwd_barcode: str
    fwd_primer: str
    rev_barcode: str
    rev_primer: str



class ReadSimulator:
    """Generates synthetic nanopore reads with controlled error rates"""

    def __init__(self,
                 noise_length: int = 100,
                 error_rate: float = 0.1,
                 homopolymer_error_rate: float = 0.2,
                 truncation_max: int = 20,
                 truncation_prob: float = 0.3,
                 adapter_seq: str = None):
        self.noise_length = noise_length
        self.error_rate = error_rate
        self.homopolymer_error_rate = homopolymer_error_rate
        self.truncation_max = truncation_max
        self.truncation_prob = truncation_prob
        self.adapter_seq = adapter_seq

        # Precompute some common values
        self.bases = ['A', 'C', 'G', 'T']
        self.error_types = ['sub', 'ins', 'del']

        # Create substitution lookup for each base
        self.sub_bases = {
            'A': ['C', 'G', 'T'],
            'C': ['A', 'G', 'T'],
            'G': ['A', 'C', 'T'],
            'T': ['A', 'C', 'G']
        }

    def weighted_truncation_length(self) -> int:
        return random.randint(0, 60)

    def truncate_sequence(self, sequence: str) -> str:
        """Apply truncation to one end of the sequence"""
        if random.random() > self.truncation_prob:
            return sequence

        trunc_len = self.weighted_truncation_length()
        if trunc_len == 0:
            return sequence

        # Randomly choose which end to truncate - bias towards final end
        if random.random() < 0.05:
            # Truncate start
            return sequence[trunc_len:]
        else:
            # Truncate end
            return sequence[:-trunc_len]

    def generate_noise(self, length: int) -> str:
        """Generate random DNA sequence of specified length"""
        return ''.join(random.choice(self.bases) for _ in range(length))

    def is_homopolymer(self, bases: list, pos: int) -> bool:
        """Check if position is in a homopolymer region (2+ identical bases)"""
        seq_len = len(bases)
        if seq_len < 2:  # Too short for homopolymers
            return False

        if pos < 0 or pos >= seq_len:  # Out of bounds check
            return False

        # Check previous base if not at start
        prev_match = pos > 0 and bases[pos] == bases[pos - 1]

        # Check next base if not at end
        next_match = pos < seq_len - 1 and bases[pos] == bases[pos + 1]

        return prev_match or next_match

    def introduce_errors(self, sequence: str, protect_start: int, protect_end: int) -> str:
        """Add random mutations, insertions, and deletions to sequence"""
        bases = list(sequence)  # Work with list instead of string
        seq_len = len(bases)

        # Validate protection bounds
        protect_start = max(0, min(protect_start, seq_len))
        protect_end = max(0, min(protect_end, seq_len))

        i = 0
        while i < seq_len:
            # Skip protected region
            if protect_start <= i < protect_end:
                i += 1
                continue

            # Check for homopolymer using list instead of string
            error_prob = (self.homopolymer_error_rate
                          if self.is_homopolymer(bases, i)
                          else self.error_rate)

            if random.random() < error_prob:
                error_type = random.choice(self.error_types)

                if error_type == 'sub' and i < seq_len:
                    # Use precomputed substitution bases
                    bases[i] = random.choice(self.sub_bases[bases[i]])
                elif error_type == 'ins':
                    bases.insert(i, random.choice(self.bases))
                    seq_len += 1
                elif error_type == 'del' and i < seq_len:
                    bases.pop(i)
                    seq_len -= 1
                    i -= 1
            i += 1

        return ''.join(bases)

    def construct_read(self, specimen: SpecimenInfo) -> str:
        """
        Construct a complete read for a specimen
        Returns (sequence, golden_string)
        """
        # Create golden string with just barcodes and primers
        golden = (f"{START_DELIMITER}"
                  f"{specimen.fwd_barcode}"
                  f"{specimen.fwd_primer}"
                  f"{specimen.rev_primer}"
                  f"{specimen.rev_barcode}"
                  f"{END_DELIMITER}")

        # Build parts of sequence separately for better memory efficiency
        parts = [
            self.adapter_seq,
            specimen.fwd_barcode,
            specimen.fwd_primer,
            self.generate_noise(self.noise_length),
            golden,
            self.generate_noise(self.noise_length),
            reverse_complement(specimen.rev_primer),
            reverse_complement(specimen.rev_barcode),
            reverse_complement(self.adapter_seq),
        ]

        # Join parts into full sequence
        sequence = ''.join(parts)

        # Add adapter and handle orientation
        is_reverse = random.choice([True, False])
        if is_reverse:
            sequence = reverse_complement(sequence)

        sequence = self.introduce_errors(sequence, ERROR_RANGE, len(sequence) - ERROR_RANGE)

        # Apply truncation
        sequence = self.truncate_sequence(sequence)

        return sequence


def read_index_file(filename: Path) -> List[SpecimenInfo]:
    """Read specimen information from Index.txt file"""
    specimens = []

    with open(filename) as f:
        reader = csv.reader(f, delimiter='\t')

        # Skip header if present
        first_row = next(reader)
        if not all(x.upper() in 'ACGT' for x in first_row[1]):  # Check if first row is header
            first_row = next(reader)

        # Process first row
        specimens.append(SpecimenInfo(
            id=first_row[0],
            fwd_barcode=first_row[1].upper(),
            fwd_primer=first_row[2].upper(),
            rev_barcode=first_row[3].upper(),
            rev_primer=first_row[4].upper()
        ))

        # Process remaining rows
        for row in reader:
            specimens.append(SpecimenInfo(
                id=row[0],
                fwd_barcode=row[1].upper(),
                fwd_primer=row[2].upper(),
                rev_barcode=row[3].upper(),
                rev_primer=row[4].upper()
            ))

    return specimens

def load_adapter_sequence(filename: str) -> str:
    """Load ONT adapter sequence from file"""
    with open(filename) as f:
        return f.read().strip()


def main():
    parser = argparse.ArgumentParser(
        description='Generate synthetic nanopore reads for barcode testing')

    parser.add_argument('index_file', type=Path,
                        help='Index.txt file containing specimen information')
    parser.add_argument('output_file', type=Path,
                        help='Output FASTQ file for simulated reads')
    parser.add_argument('--reads-per-specimen', type=int, default=200,
                        help='Number of reads to generate per specimen')
    parser.add_argument('--noise-length', type=int, default=100,
                        help='Length of random noise sections')
    parser.add_argument('--error-rate', type=float, default=0.10,
                        help='Base error rate for mutations')
    parser.add_argument('--homopolymer-error-rate', type=float, default=0.10,
                        help='Error rate for homopolymer regions')
    parser.add_argument('--truncation-max', type=int, default=20,
                        help='Maximum bases to truncate from sequence end')
    parser.add_argument('--truncation-prob', type=float, default=0.5,
                        help='Probability of truncating sequence')
    parser.add_argument('--adapter-file', type=str,
                        help='File containing ONT adapter sequence')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility')
    parser.add_argument('--batch-size', type=int, default=1000,
                        help='Number of reads to write in each batch')

    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Set random seed if provided
    if args.seed is not None:
        random.seed(args.seed)

    # Load adapter sequence
    adapter_seq = None
    if args.adapter_file:
        adapter_seq = load_adapter_sequence(args.adapter_file)
    else:
        adapter_seq = "CCTGTACTTCGTTCAGTTACGTATTGCT"  # Default ONT adapter
        logging.warning("Using default ONT adapter sequence. Consider providing adapter file.")

    # Initialize simulator
    simulator = ReadSimulator(
        noise_length=args.noise_length,
        error_rate=args.error_rate,
        homopolymer_error_rate=args.homopolymer_error_rate,
        truncation_max=args.truncation_max,
        truncation_prob=args.truncation_prob,
        adapter_seq=adapter_seq
    )

    # Read specimen information
    specimens = read_index_file(args.index_file)
    logging.info(f"Loaded {len(specimens)} specimens from index file")

    total_reads = len(specimens) * args.reads_per_specimen
    read_count = 0
    batch = []
    progress_interval = max(1, total_reads // 100)  # Update progress every 1%

    with open(args.output_file, 'w') as f:
        print("Generating reads: ", end='', flush=True)

        for specimen in specimens:
            for i in range(args.reads_per_specimen):
                sequence = simulator.construct_read(specimen)
                read_id = f"{specimen.id}_read_{i}"

                # Build FASTQ entry
                batch.extend([
                    f"@{read_id}",
                    sequence,
                    "+",
                    "~" * len(sequence)
                ])

                read_count += 1

                # Update progress bar
                if read_count % progress_interval == 0:
                    progress = read_count / total_reads * 100
                    print(f"\rGenerating reads: {progress:3.0f}%", end='', flush=True)

                # Write batch if it's full
                if len(batch) >= args.batch_size * 4:
                    f.write('\n'.join(batch) + '\n')
                    batch = []

        # Write any remaining reads
        if batch:
            f.write('\n'.join(batch) + '\n')

        print("\rGenerating reads: 100%")  # Final progress update

    logging.info(f"Generated {read_count} reads")

if __name__ == "__main__":
    main()

