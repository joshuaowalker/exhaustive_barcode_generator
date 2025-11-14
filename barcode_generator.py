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
class DimerConfig:
    """Configuration for self-dimer checking"""
    min_overlap: int = 4
    score_threshold: float = 0.35
    continuous_weight: float = 0.7
    total_weight: float = 0.3

    def __post_init__(self):
        if not 0 <= self.score_threshold <= 1:
            raise ValueError("score_threshold must be between 0 and 1")
        if not 0 <= self.continuous_weight <= 1:
            raise ValueError("continuous_weight must be between 0 and 1")
        if not 0 <= self.total_weight <= 1:
            raise ValueError("total_weight must be between 0 and 1")
        if abs(self.continuous_weight + self.total_weight - 1.0) > 1e-6:
            raise ValueError("continuous_weight and total_weight must sum to 1")

class SelfDimerChecker:
    """Checks DNA sequences for self-dimerization potential"""

    def __init__(self, config: DimerConfig):
        self.config = config

    @staticmethod
    def _is_complementary(base1: str, base2: str) -> bool:
        """Check if two bases are complementary"""
        pairs = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C'
        }
        return pairs.get(base1) == base2

    def _score_alignment(self, matches: list) -> float:
        """Score an alignment based on match quality
        Returns a score between 0 and 1, with 1 being strongest interaction"""
        if not matches:
            return 0.0

        # Count continuous matches
        max_continuous = 0
        current_continuous = 0
        total_matches = 0

        for match in matches:
            if match:
                current_continuous += 1
                total_matches += 1
                max_continuous = max(max_continuous, current_continuous)
            else:
                current_continuous = 0

        # Weight continuous matches more heavily using configured weights
        score = (self.config.continuous_weight * (max_continuous / len(matches))) + \
                (self.config.total_weight * (total_matches / len(matches)))

        return score

    def check_self_dimers(self, sequence: str) -> tuple[bool, str]:
        """Check sequence for self-dimerization potential
        Returns (has_dimers, alignment_str)"""
        sequence = sequence.upper()
        best_score = 0
        best_alignment = ""

        # Try different alignments/overlaps
        for offset in range(-len(sequence) + self.config.min_overlap,
                            len(sequence) - self.config.min_overlap + 1):
            matches = []
            alignment_str = ""

            # Create alignment strings
            top = sequence
            bottom = sequence[::-1]  # Reverse for complementary check

            if offset < 0:
                top = " " * abs(offset) + top
            else:
                bottom = " " * offset + bottom

            # Trim to matching lengths
            max_len = max(len(top), len(bottom))
            top = top.ljust(max_len)
            bottom = bottom.ljust(max_len)

            # Check complementary bases
            for i in range(min(len(top), len(bottom))):
                if top[i] != " " and bottom[i] != " ":
                    is_match = SelfDimerChecker._is_complementary(top[i], bottom[i])
                    matches.append(is_match)
                    alignment_str += "|" if is_match else " "
                else:
                    alignment_str += " "

            # Score this alignment
            score = self._score_alignment(matches)

            # Keep track of best alignment
            if score > best_score:
                best_score = score
                best_alignment = f"5-{sequence}->\n{alignment_str}\n   <-{sequence[::-1]}-5"

        has_dimers = best_score > self.config.score_threshold
        return has_dimers, best_alignment


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
    filter_rc: bool
    current_barcodes: List[List[str]]
    exclusion_seqs: List[str]
    dimer_checker: Optional[SelfDimerChecker] = None  # Add this line




class ParallelBarcodeGenerator:
    """Manages parallel generation and validation of barcodes"""

    def __init__(self,
                 length: int,
                 min_distance: int,
                 min_gc: float,
                 max_gc: float,
                 max_homopolymer: int,
                 filter_rc: bool = False,
                 dimer_config: Optional[DimerConfig] = None,
                 exclusion_file: Optional[Path] = None,
                 initial_barcodes: Optional[Path] = None,  # New parameter
                 random_seqs: int = 1,
                 num_processes: Optional[int] = None,
                 chunk_size: int = 10000):

        self.length = length
        self.min_distance = min_distance
        self.min_gc = min_gc
        self.max_gc = max_gc
        self.max_homopolymer = max_homopolymer
        self.filter_rc = filter_rc
        self.dimer_config = dimer_config or DimerConfig()
        self.dimer_checker = SelfDimerChecker(self.dimer_config)
        self.num_processes = num_processes or cpu_count()
        self.chunk_size = chunk_size
        self.bases = ['A', 'C', 'G', 'T']
        self.total_barcodes = 4 ** length
        self.exclusion_seqs = self._load_exclusion_seqs(exclusion_file, random_seqs)

        logging.info(f"Using {self.num_processes} processes to check {self.total_barcodes:,} sequences")
        logging.info(f"Using chunk size of {self.chunk_size:,} sequences per work item")

        logging.info(f"Using {self.num_processes} processes to check {self.total_barcodes:,} sequences")
        logging.info(f"Using chunk size of {self.chunk_size:,} sequences per work item")

        # Load initial barcodes if provided
        self.initial_barcodes = []
        if initial_barcodes:
            self.initial_barcodes = load_initial_barcodes(initial_barcodes, length)
            logging.info(f"Loaded {len(self.initial_barcodes)} initial barcodes")

            # Validate initial barcodes against constraints
            valid_barcodes = []
            for barcode in self.initial_barcodes:
                barcode_str = ''.join(barcode)
                gc = self.gc_cont(barcode)
                has_issues = False

                if not (min_gc <= gc <= max_gc):
                    logging.warning(f"Initial barcode {barcode_str} has invalid GC content: {gc:.2%}")
                    has_issues = True

                if self.has_homopolymer(barcode, max_homopolymer):
                    logging.warning(f"Initial barcode {barcode_str} has homopolymer longer than {max_homopolymer}")
                    has_issues = True

                # Check distance with previously validated barcodes
                is_compatible = True
                for valid_barcode in valid_barcodes:
                    if not self.is_compatible(barcode, [valid_barcode],
                                              self.exclusion_seqs, min_distance, filter_rc)[0]:
                        logging.warning(
                            f"Initial barcode {barcode_str} violates minimum distance requirement with {''.join(valid_barcode)}")
                        has_issues = True
                        is_compatible = False
                        break

                if not has_issues and is_compatible:
                    valid_barcodes.append(barcode)
                else:
                    logging.warning(f"Skipping initial barcode {barcode_str} due to constraint violations")

            if len(valid_barcodes) < len(self.initial_barcodes):
                logging.warning(
                    f"Only {len(valid_barcodes)} of {len(self.initial_barcodes)} initial barcodes met all constraints")

            self.initial_barcodes = valid_barcodes

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
                      min_distance: int,
                      filter_rc: bool) -> Tuple[bool, int]:
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
            elif filter_rc and edlib.align(new_str_rc, existing_str, task='distance')['editDistance'] < min_distance:
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

        idx = work_item.start_idx
        while idx < min(work_item.start_idx + work_item.chunk_size, 4 ** work_item.length):
            barcode = ParallelBarcodeGenerator.index_to_barcode(idx, work_item.length)

            # Check GC content first as it's fast
            if not work_item.min_gc <= ParallelBarcodeGenerator.gc_cont(barcode) <= work_item.max_gc:
                idx += 1
                continue

            # Check homopolymers with skip distance
            has_homo, skip = ParallelBarcodeGenerator.has_homopolymer_with_skip(
                barcode, work_item.max_homopolymer)
            if has_homo:
                idx += max(1, skip)  # Always advance at least 1
                continue

            # Check compatibility with current set
            compatible, conflict_idx = ParallelBarcodeGenerator.is_compatible(
                barcode, local_barcodes, work_item.exclusion_seqs,
                work_item.min_distance, work_item.filter_rc)

            if compatible:
                # Post checks - these are potentially slower
                sequence = "".join(barcode)
                if work_item.dimer_checker:
                    has_dimers, _ = work_item.dimer_checker.check_self_dimers(sequence)
                    if has_dimers:
                        idx += 1
                        continue

                if ParallelBarcodeGenerator.has_problematic_repeats(sequence):
                    idx += 1
                    continue

                candidate_barcodes.append(barcode)
            elif conflict_idx >= 0:
                # Move conflicting barcode to front of list for faster future checks
                conflicting = local_barcodes.pop(conflict_idx)
                local_barcodes.insert(0, conflicting)

            idx += 1

        return candidate_barcodes

    @staticmethod
    def calculate_skip_distance(barcode: List[str], position: int, max_homopolymer: int) -> int:
        """
        Calculate how many sequences to skip based on a homopolymer violation.
        Returns the number of sequences to skip to clear the homopolymer.

        Args:
            barcode: The current barcode sequence
            position: Position where homopolymer violation was detected
            max_homopolymer: Maximum allowed homopolymer length
        """
        # Count the homopolymer length up to this position
        homo_len = 1
        base = barcode[position]
        for i in range(position - 1, -1, -1):
            if barcode[i] != base:
                break
            homo_len += 1

        # Calculate how many positions we need to skip to clear this homopolymer
        excess = homo_len - max_homopolymer
        if excess <= 0:
            return 0

        # Calculate skip distance based on position
        # Each position represents 4^n sequences where n is the remaining length
        remaining_length = len(barcode) - position - 1
        skip_distance = 4 ** remaining_length

        return skip_distance

    @staticmethod
    def has_homopolymer_with_skip(barcode: List[str], max_homopolymer_length: int) -> Tuple[bool, int]:
        """
        Modified homopolymer check that returns both the result and how many sequences to skip.
        Returns (has_homopolymer, skip_distance)
        """
        count = 1
        for i in range(1, len(barcode)):
            if barcode[i] == barcode[i - 1]:
                count += 1
                if count > max_homopolymer_length:
                    # Calculate skip distance when violation is found
                    skip = ParallelBarcodeGenerator.calculate_skip_distance(
                        barcode, i, max_homopolymer_length)
                    return True, skip
            else:
                count = 1
        return False, 0

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

    @staticmethod
    def detect_adjacent_repeats(sequence: str, min_repeat_len: int = 2, max_repeat_len: int = 4) -> List[
        Tuple[str, int]]:
        """
        Find adjacent repeating patterns of length min_repeat_len to max_repeat_len.
        Returns list of (pattern, count) tuples for found repeats.

        Example:
            "TCTGTGTCTG" with min_repeat_len=2, max_repeat_len=3 would find:
            [("TG", 2)] for the adjacent TGTG
        """
        repeats = []

        # Check each possible repeat length
        for length in range(min_repeat_len, max_repeat_len + 1):
            i = 0
            while i <= len(sequence) - length:
                pattern = sequence[i:i + length]
                # Count consecutive occurrences
                count = 1
                pos = i + length
                while pos <= len(sequence) - length and sequence[pos:pos + length] == pattern:
                    count += 1
                    pos += length

                if count > 1:
                    repeats.append((pattern, count))
                    i = pos  # Skip past this repeat
                else:
                    i += 1

        return repeats

    @staticmethod
    def has_problematic_repeats(sequence: str,
                                min_repeat_len: int = 2,
                                max_repeat_len: int = 4,
                                max_repeat_count: int = 2) -> bool:
        """
        Check if sequence has problematic adjacent repeats.
        Returns True if problematic repeats found, False otherwise.
        """
        repeats = ParallelBarcodeGenerator.detect_adjacent_repeats(sequence, min_repeat_len, max_repeat_len)

        for pattern, count in repeats:
            if count >= max_repeat_count:
                return True

        return False

    def generate_barcodes(self, target_size: int) -> Tuple[List[List[str]], int]:
        """Generate barcodes using producer/consumer pattern"""
        # Start with initial barcodes if provided
        accepted_barcodes = self.initial_barcodes.copy()
        sequences_checked = 0
        next_index = 0

        # Adjust target size to account for initial barcodes
        remaining_target = max(0, target_size - len(accepted_barcodes))

        if len(accepted_barcodes) >= target_size:
            logging.info(f"Initial barcode set already meets target size ({len(accepted_barcodes)} >= {target_size})")
            return accepted_barcodes, sequences_checked

        with Pool(processes=self.num_processes) as pool:
            pbar = tqdm(total=self.total_barcodes, desc="Testing sequences")

            while sequences_checked < self.total_barcodes and len(accepted_barcodes) < target_size:
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
                        filter_rc=self.filter_rc,
                        current_barcodes=accepted_barcodes.copy(),
                        exclusion_seqs=self.exclusion_seqs,
                        dimer_checker=self.dimer_checker  # Add this line
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
                                                           self.min_distance, self.filter_rc)
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

def load_initial_barcodes(file_path: Path, length: int) -> List[List[str]]:
    """Load pre-specified barcodes from file and validate them"""
    barcodes = []

    with open(file_path) as f:
        for line_num, line in enumerate(f, 1):
            # Skip empty lines and comments
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Extract barcode sequence (ignore any additional fields)
            barcode = line.split()[0].upper()

            # Validate barcode
            if len(barcode) != length:
                raise ValueError(f"Barcode on line {line_num} has incorrect length: {barcode}")
            if not set(barcode).issubset({'A', 'C', 'G', 'T'}):
                raise ValueError(f"Barcode on line {line_num} contains invalid characters: {barcode}")

            # Convert to list format
            barcodes.append(list(barcode))

    return barcodes



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
    parser.add_argument(
        '--filter-rc',
        action='store_true',
        help='Filter barcodes with reverse complement conflicts'
    )
    parser.add_argument(
        '--initial-barcodes',
        type=Path,
        help='File containing pre-specified barcodes (one per line)'
    )

    parser.add_argument(
        '--dimer-min-overlap',
        type=int,
        default=4,
        help='Minimum overlap length to consider for self-dimers (default: 4)'
    )
    parser.add_argument(
        '--dimer-score-threshold',
        type=float,
        default=0.35,
        help='Score threshold above which to consider a sequence as having self-dimers (0-1, default: 0.35)'
    )
    parser.add_argument(
        '--dimer-continuous-weight',
        type=float,
        default=0.7,
        help='Weight given to continuous matches in dimer scoring (0-1, default: 0.7)'
    )
    parser.add_argument(
        '--dimer-total-weight',
        type=float,
        default=0.3,
        help='Weight given to total matches in dimer scoring (0-1, default: 0.3)'
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
        f.write(f'Found {len(barcodes)} barcodes:\n\n')

        # Write header
        f.write('Sequence\tGC Content\tHomopolymer %\n')

        # Write barcodes with metrics
        for barcode_str, gc, homo_pct in barcode_data:
            f.write(f'{barcode_str}\t{gc:.2%}\t{homo_pct:.2%}\n')


def main() -> None:
    args = parse_arguments()
    setup_logging(args.debug)

    # Create dimer configuration
    dimer_config = DimerConfig(
        min_overlap=args.dimer_min_overlap,
        score_threshold=args.dimer_score_threshold,
        continuous_weight=args.dimer_continuous_weight,
        total_weight=args.dimer_total_weight
    )

    # Initialize generator
    generator = ParallelBarcodeGenerator(
        length=args.length,
        min_distance=args.min_distance,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        max_homopolymer=args.max_homopolymer,
        filter_rc=args.filter_rc,
        dimer_config=dimer_config,
        exclusion_file=args.exclusion_file,
        initial_barcodes=args.initial_barcodes,  # Add this line
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


    # Write output with sorting
    write_output(
        barcodes,
        args.output,
        args.length,
        args.min_distance,
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
    print("\nGenerated barcodes with metrics:")
    print("Sequence\tGC Content\tHomopolymer %")
    for barcode_str, gc, homo_pct in barcode_data:
        print(f"{barcode_str}\t{gc:.2%}\t{homo_pct:.2%}")

    logging.info('Complete')

if __name__ == '__main__':
    main()


