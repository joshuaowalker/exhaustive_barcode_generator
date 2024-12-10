#! /usr/bin/env python3

'''PARALLEL BARCODE GENERATOR
   Systematically generates DNA barcodes using producer/consumer pattern'''

import sys
import time
import re
import logging
import argparse
from datetime import date
from pathlib import Path
from typing import List, Tuple, Optional, Iterator, Dict, Set
from Bio import SeqIO
from Bio.Seq import reverse_complement
import edlib
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from dataclasses import dataclass
from contextlib import contextmanager
from functools import partial


@dataclass
class WorkItem:
    """Parameters for worker processes"""
    start_idx: int
    chunk_size: int
    length: int
    min_distance: int
    min_gc: float
    max_gc: float
    max_homopolymer: int
    current_barcodes: List[List[str]]
    exclusion_seqs: List[str]


class ParallelBarcodeGenerator:
    """Manages parallel generation and validation of barcodes"""

    def __init__(self,
                 length: int,
                 min_distance: int,
                 min_gc: float,
                 max_gc: float,
                 max_homopolymer: int,
                 exclusion_file: Optional[Path] = None,
                 random_seqs: int = 1,
                 num_processes: Optional[int] = None,
                 chunk_size: int = 10000):
        self.length = length
        self.min_distance = min_distance
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.max_homopolymer = max_homopolymer
        self.num_processes = num_processes or cpu_count()
        self.chunk_size = chunk_size
        self.bases = ['A', 'C', 'G', 'T']
        self.total_barcodes = 4 ** length
        self.exclusion_seqs = self._load_exclusion_seqs(exclusion_file, random_seqs)

        logging.info(f"Using {self.num_processes} processes to check {self.total_barcodes:,} sequences")
        logging.info(f"Using chunk size of {self.chunk_size:,} sequences per work item")

        logging.info(f"Using {self.num_processes} processes to check {self.total_barcodes:,} sequences")
        logging.info(f"Using chunk size of {self.chunk_size:,} sequences per work item")

    @staticmethod
    def index_to_barcode(idx: int, length: int) -> List[str]:
        """Convert an index to a barcode sequence"""
        bases = ['A', 'C', 'G', 'T']
        barcode = []
        for j in range(length - 1, -1, -1):
            base_idx = (idx >> (j * 2)) & 0b11
            barcode.append(bases[base_idx])
        return barcode

    @staticmethod
    def gc_cont(barcode: List[str]) -> float:
        """Calculate GC content of a barcode"""
        gc = sum(1 for base in barcode if base in ['C', 'G'])
        return gc / len(barcode)

    @staticmethod
    def has_homopolymer(barcode: List[str], max_homopolymer_length: int) -> bool:
        """Check if barcode has homopolymer runs longer than specified length"""
        count = 1
        for i in range(1, len(barcode)):
            if barcode[i] == barcode[i - 1]:
                count += 1
                if count > max_homopolymer_length:
                    return True
            else:
                count = 1
        return False

    @staticmethod
    def is_compatible(new_barcode: List[str],
                      accepted_barcodes: List[List[str]],
                      exclusion_seqs: List[str],
                      min_distance: int) -> Tuple[bool, int]:
        """
        Check if a new barcode meets the distance requirement with all accepted barcodes
        Returns (is_compatible, index_of_conflict) where index_of_conflict is -1 if compatible
        """
        new_str = ''.join(new_barcode)
        new_str_rc = reverse_complement(new_str)

        # Check against accepted barcodes
        for i, existing in enumerate(accepted_barcodes):
            existing_str = ''.join(existing)
            if edlib.align(new_str, existing_str, task='distance')['editDistance'] < min_distance:
                return False, i
            elif edlib.align(new_str_rc, existing_str, task='distance')['editDistance'] < min_distance:
                return False, i

        # Check against exclusion sequences
        for seq in exclusion_seqs:
            if edlib.align(new_str, seq, task='distance')['editDistance'] < min_distance:
                return False, -1
            elif edlib.align(new_str_rc, seq, task='distance')['editDistance'] < min_distance:
                return False, -1

        # Check against self (reverse complement)
        if edlib.align(new_str, new_str_rc, task='distance')['editDistance'] < min_distance:
            return False, -1

        return True, -1

    @staticmethod
    def worker_process(work_item: WorkItem) -> List[List[str]]:
        """Worker process for barcode generation and validation"""
        candidate_barcodes = []
        local_barcodes = work_item.current_barcodes.copy()

        for idx in range(work_item.start_idx, work_item.start_idx + work_item.chunk_size):
            if idx >= 4 ** work_item.length:
                break

            barcode = ParallelBarcodeGenerator.index_to_barcode(idx, work_item.length)

            # Early filtering
            if not (work_item.min_gc <= ParallelBarcodeGenerator.gc_cont(barcode) <= work_item.max_gc and
                    not ParallelBarcodeGenerator.has_homopolymer(barcode, work_item.max_homopolymer)):
                continue

            # Check compatibility with current set
            compatible, conflict_idx = ParallelBarcodeGenerator.is_compatible(
                barcode, local_barcodes, work_item.exclusion_seqs, work_item.min_distance)

            if compatible:
                candidate_barcodes.append(barcode)
            elif conflict_idx >= 0:
                # Move conflicting barcode to front of list for faster future checks
                conflicting = local_barcodes.pop(conflict_idx)
                local_barcodes.insert(0, conflicting)

        return candidate_barcodes

    @staticmethod
    def homopolymer_percentage(barcode: List[str]) -> float:
        """Calculate the percentage of bases that are part of homopolymer runs"""
        if not barcode:
            return 0.0

        homopolymer_positions = 0
        current_run = 1

        # Count positions in homopolymer runs
        for i in range(1, len(barcode)):
            if barcode[i] == barcode[i - 1]:
                current_run += 1
                # Count both positions in a run of 2, and all positions in longer runs
                if current_run == 2:
                    homopolymer_positions += 2
                elif current_run > 2:
                    homopolymer_positions += 1
            else:
                current_run = 1

        return homopolymer_positions / len(barcode)

    def _generate_random_sequences(self, n: int = 1) -> List[str]:
        """Generate n random sequences of barcode length"""
        import random
        seqs = []
        bases = ['A', 'C', 'G', 'T']
        for _ in range(n):
            seq = ''.join(random.choice(bases) for _ in range(self.length))
            seqs.append(seq)
        return seqs

    def _load_exclusion_seqs(self, filename: Optional[Path], random_seqs: int = 1) -> List[str]:
        """Load and process exclusion sequences from FASTA file"""
        exclusion_seqs = []

        # Add random sequences first
        if random_seqs > 0:
            rand_seqs = self._generate_random_sequences(random_seqs)
            exclusion_seqs.extend(rand_seqs)
            logging.info(f"Added {len(rand_seqs)} random sequences to exclusion set")

        if filename is None:
            return exclusion_seqs

        iupac_codes = set('RYSWKMBDHVN')

        with open(filename) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper()

                # Check for IUPAC codes
                if any(base in iupac_codes for base in seq):
                    logging.warning(f"IUPAC ambiguity codes found in sequence {record.id}")
                    continue

                # Add sequence and its reverse complement
                processed_seqs = []

                # Handle sequences shorter than barcode length
                if len(seq) < self.length:
                    padded = seq + 'A' * (self.length - len(seq))
                    processed_seqs.append(padded)
                # Handle sequences longer than barcode length
                elif len(seq) > self.length:
                    for i in range(len(seq) - self.length + 1):
                        subseq = seq[i:i + self.length]
                        processed_seqs.append(subseq)
                else:
                    processed_seqs.append(seq)

                # Add sequences and their reverse complements
                for seq in processed_seqs:
                    exclusion_seqs.append(seq)

        logging.info(f"Loaded {len(exclusion_seqs)} exclusion sequences")
        return exclusion_seqs

    def generate_barcodes(self, target_size: int) -> Tuple[List[List[str]], int]:
        """Generate barcodes using producer/consumer pattern"""
        accepted_barcodes = []
        sequences_checked = 0
        next_index = 0

        with Pool(processes=self.num_processes) as pool:
            pbar = tqdm(total=self.total_barcodes, desc="Testing sequences")

            while sequences_checked < self.total_barcodes and len(accepted_barcodes) < target_size:
                # Create work items
                work_items = []
                for _ in range(self.num_processes):
                    if next_index >= self.total_barcodes:
                        break

                    work_items.append(WorkItem(
                        start_idx=next_index,
                        chunk_size=self.chunk_size,
                        length=self.length,
                        min_distance=self.min_distance,
                        min_gc=self.min_gc,
                        max_gc=self.max_gc,
                        max_homopolymer=self.max_homopolymer,
                        current_barcodes=accepted_barcodes.copy(),
                        exclusion_seqs=self.exclusion_seqs
                    ))
                    next_index += self.chunk_size

                if not work_items:
                    break

                    # Process work items in parallel
                for candidate_barcodes in pool.imap_unordered(self.worker_process, work_items):
                    # Validate results sequentially to handle race conditions
                    for barcode in candidate_barcodes:
                        if len(accepted_barcodes) >= target_size:
                            break
                        compatible, _ = self.is_compatible(barcode, accepted_barcodes, self.exclusion_seqs,
                                                           self.min_distance)
                        if compatible:
                            accepted_barcodes.append(barcode)

                    # Update progress
                    sequences_checked += self.chunk_size
                    pbar.set_postfix({'found': len(accepted_barcodes)}, refresh=True)
                    pbar.update(self.chunk_size)

                if len(accepted_barcodes) >= target_size:
                    break

            pbar.close()

        return accepted_barcodes, sequences_checked


def setup_logging(debug: bool = False) -> None:
    """Configure logging settings"""
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Generate DNA barcodes using parallel systematic enumeration',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '--length',
        type=int,
        default=13,
        help='Length of barcodes to generate'
    )
    parser.add_argument(
        '--number',
        type=int,
        default=96,
        help='Number of barcodes to generate'
    )
    parser.add_argument(
        '--min-distance',
        type=int,
        default=7,
        help='Minimum Levenshtein distance between barcodes'
    )
    parser.add_argument(
        '--min-gc',
        type=float,
        default=0.4,
        help='Minimum GC content (0-1)'
    )
    parser.add_argument(
        '--max-gc',
        type=float,
        default=0.6,
        help='Maximum GC content (0-1)'
    )
    parser.add_argument(
        '--max-homopolymer',
        type=int,
        default=2,
        help='Maximum allowed homopolymer length'
    )
    parser.add_argument(
        '--processes',
        type=int,
        default=None,
        help='Number of worker processes (default: CPU count)'
    )
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=10000,
        help='Number of sequences to process per work item'
    )
    parser.add_argument(
        '--output',
        type=Path,
        help='Output file path (default: barcode_<date>.txt)'
    )
    parser.add_argument(
        '--sort',
        choices=['homopolymer', 'gc_content', 'lexicographic'],
        default='homopolymer',
        help='Sort output by criterion (homopolymer: ascending order of homopolymer percentage)'
    )
    parser.add_argument(
        '--random-seqs',
        type=int,
        default=1,
        help='Number of random sequences to add to exclusion set'
    )
    parser.add_argument(
        '--exclusion-file',
        type=Path,
        help='FASTA file containing sequences to exclude (e.g., primers, adapters)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging'
    )

    args = parser.parse_args()

    # Validate arguments
    if args.min_gc < 0 or args.max_gc > 1:
        parser.error("GC content must be between 0 and 1")
    if args.min_gc > args.max_gc:
        parser.error("min-gc must be less than max-gc")
    if args.length < 1:
        parser.error("Length must be positive")
    if args.number < 1:
        parser.error("Number must be positive")
    if args.min_distance < 1:
        parser.error("Minimum distance must be positive")
    if args.max_homopolymer < 1:
        parser.error("Maximum homopolymer length must be positive")
    if args.processes is not None and args.processes < 1:
        parser.error("Number of processes must be positive")
    if args.chunk_size < 1:
        parser.error("Chunk size must be positive")

    # Set default output path if not provided
    if args.output is None:
        args.output = Path(f'barcode_{date.today().isoformat()}.txt')

    return args

def write_output(
        barcodes: List[List[str]],
        output_path: Path,
        length: int,
        min_distance: int,
        runtime: float,
        sequences_checked: int,
        sort_by: Optional[str] = None
) -> None:
    """Write barcodes and statistics to output file"""

    # Prepare barcode data with metrics
    barcode_data = []
    for barcode in barcodes:
        barcode_str = ''.join(barcode)
        gc = ParallelBarcodeGenerator.gc_cont(barcode)
        homo_pct = ParallelBarcodeGenerator.homopolymer_percentage(barcode)
        barcode_data.append((barcode_str, gc, homo_pct))

    # Sort if requested
    if sort_by == 'gc_content':
        barcode_data.sort(key=lambda x: x[1])
    elif sort_by == 'lexicographic':
        barcode_data.sort(key=lambda x: x[0])
    elif sort_by == 'homopolymer':
        barcode_data.sort(key=lambda x: x[2])

    with open(output_path, 'w') as f:
        f.write(f'Parallel Barcode Generator Results\n')
        f.write(f'Generated on: {date.today().isoformat()}\n\n')
        f.write(f'Parameters:\n')
        f.write(f'- Barcode length: {length}\n')
        f.write(f'- Minimum Levenshtein distance: {min_distance}\n')
        f.write(f'- Sequences checked: {sequences_checked:,}\n')
        f.write(f'- Runtime: {runtime:.2f} seconds\n\n')
        f.write(f'Found {len(barcodes)} barcodes:\n\n')

        # Write header
        f.write('Sequence\tGC Content\tHomopolymer %\n')

        # Write barcodes with metrics
        for barcode_str, gc, homo_pct in barcode_data:
            f.write(f'{barcode_str}\t{gc:.2%}\t{homo_pct:.2%}\n')



def main() -> None:
    """Main function"""
    args = parse_arguments()
    setup_logging(args.debug)

    start_time = time.process_time()

    # Initialize generator
    generator = ParallelBarcodeGenerator(
        length=args.length,
        min_distance=args.min_distance,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        max_homopolymer=args.max_homopolymer,
        exclusion_file=args.exclusion_file,
        random_seqs=args.random_seqs,
        num_processes=args.processes,
        chunk_size=args.chunk_size
    )

    # Generate barcodes
    try:
        barcodes, sequences_checked = generator.generate_barcodes(args.number)
    except KeyboardInterrupt:
        logging.info("\nInterrupted by user. Saving results...")
        barcodes, sequences_checked = [], 0

    end_time = time.process_time()
    runtime = end_time - start_time

    # Write output with sorting
    write_output(
        barcodes,
        args.output,
        args.length,
        args.min_distance,
        runtime,
        sequences_checked,
        args.sort
    )

    # Prepare barcode data with metrics
    barcode_data = []
    for barcode in barcodes:
        barcode_str = ''.join(barcode)
        gc = ParallelBarcodeGenerator.gc_cont(barcode)
        homo_pct = ParallelBarcodeGenerator.homopolymer_percentage(barcode)
        barcode_data.append((barcode_str, gc, homo_pct))

    # Sort if requested
    if args.sort == 'homopolymer':
        barcode_data.sort(key=lambda x: x[2])  # Sort by homopolymer percentage

    # Print summary to console
    print("\nRESULTS\n")
    print(f"Sequences checked: {sequences_checked:,}")
    print(f"Valid barcodes found: {len(barcodes)}")
    print(f"Runtime: {runtime:.2f} seconds")
    print("\nGenerated barcodes with metrics:")
    print("Sequence\tGC Content\tHomopolymer %")
    for barcode_str, gc, homo_pct in barcode_data:
        print(f"{barcode_str}\t{gc:.2%}\t{homo_pct:.2%}")

    logging.info('Complete')

if __name__ == '__main__':
    main()


