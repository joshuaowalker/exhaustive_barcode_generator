#!/usr/bin/env python3

"""
Demultiplexing Validator for testing barcode demultiplexing accuracy
"""

import argparse
import csv
import re
import logging
import statistics
from pathlib import Path
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from typing import Dict, Set, List, Tuple, Optional
from Bio import SeqIO
from Bio.Seq import reverse_complement

# Golden string delimiters from simulator
START_DELIMITER = "AAAAAAAAAAGGGGGGGGGG"
END_DELIMITER =   "GGGGGGGGGGAAAAAAAAAA"


@dataclass
class ValidationStats:
    """Container for validation statistics"""
    total_reads: int = 0
    correct_assignments: int = 0
    incorrect_assignments: int = 0
    unknown_reads: int = 0
    ambiguous_reads: int = 0
    failed_golden_extraction: int = 0

    def accuracy(self) -> float:
        """Calculate percentage of correct assignments"""
        if self.total_reads == 0:
            return 0.0
        return self.correct_assignments / self.total_reads * 100

    def contamination_rate(self) -> float:
        """Calculate percentage of incorrect assignments"""
        if self.total_reads == 0:
            return 0.0
        return self.incorrect_assignments / self.total_reads * 100

    def unassigned_rate(self) -> float:
        """Calculate percentage of unassigned reads (unknown + ambiguous)"""
        if self.total_reads == 0:
            return 0.0
        return (self.unknown_reads + self.ambiguous_reads) / self.total_reads * 100


@dataclass
class BarcodeStats:
    """Statistics for individual barcodes"""
    total_reads: int = 0
    correct_reads: int = 0
    incorrect_reads: int = 0
    unassigned_reads: int = 0
    edit_distances: Counter = field(default_factory=Counter)

    def accuracy(self) -> float:
        """Calculate percentage of correct reads"""
        if self.total_reads == 0:
            return 0.0
        return self.correct_reads / self.total_reads * 100

    def contamination_rate(self) -> float:
        """Calculate percentage of incorrect reads"""
        if self.total_reads == 0:
            return 0.0
        return self.incorrect_reads / self.total_reads * 100

    def unassigned_rate(self) -> float:
        """Calculate percentage of unassigned reads"""
        if self.total_reads == 0:
            return 0.0
        return self.unassigned_reads / self.total_reads * 100

    def mean_edit_distance(self) -> float:
        """Calculate mean edit distance"""
        if not self.edit_distances:
            return 0.0
        total = sum(dist * count for dist, count in self.edit_distances.items())
        count = sum(self.edit_distances.values())
        return total / count



@dataclass
class SpecimenInfo:
    """Container for specimen information"""
    id: str
    fwd_barcode: str
    fwd_primer: str
    rev_barcode: str
    rev_primer: str


class DemuxValidator:
    """Validates demultiplexing results against golden strings"""

    def __init__(self):
        self.global_stats = ValidationStats()
        self.specimen_stats: Dict[str, ValidationStats] = defaultdict(ValidationStats)
        self.barcode_stats: Dict[str, BarcodeStats] = defaultdict(BarcodeStats)
        self.specimens: Dict[str, SpecimenInfo] = {}

        # Initialize component lengths
        self.fwd_barcode_len = None
        self.rev_barcode_len = None
        self.fwd_primer_len = None
        self.rev_primer_len = None

    def read_index_file(self, index_file: Path):
        """Read specimen information from Index.txt"""
        with open(index_file) as f:
            reader = csv.reader(f, delimiter='\t')

            # Skip header if present
            first_row = next(reader)
            if not all(x.upper() in 'ACGT' for x in first_row[1]):
                first_row = next(reader)

            # Process first row to get barcode lengths
            spec = SpecimenInfo(
                id=first_row[0],
                fwd_barcode=first_row[1].upper(),
                fwd_primer=first_row[2].upper(),
                rev_barcode=first_row[3].upper(),
                rev_primer=first_row[4].upper()
            )
            self.fwd_barcode_len = len(spec.fwd_barcode)
            self.rev_barcode_len = len(spec.rev_barcode)
            self.fwd_primer_len = len(spec.fwd_primer)
            self.rev_primer_len = len(spec.rev_primer)
            self.specimens[spec.id] = spec

            # Process remaining rows and validate consistent lengths
            for row in reader:
                spec = SpecimenInfo(
                    id=row[0],
                    fwd_barcode=row[1].upper(),
                    fwd_primer=row[2].upper(),
                    rev_barcode=row[3].upper(),
                    rev_primer=row[4].upper()
                )

                # Verify consistent barcode lengths
                if len(spec.fwd_barcode) != self.fwd_barcode_len:
                    raise ValueError(f"Inconsistent forward barcode length for {spec.id}")
                if len(spec.rev_barcode) != self.rev_barcode_len:
                    raise ValueError(f"Inconsistent reverse barcode length for {spec.id}")

                self.specimens[spec.id] = spec

                # Initialize barcode stats
                self.barcode_stats[spec.fwd_barcode]
                self.barcode_stats[spec.rev_barcode]

    def extract_golden_info(self, record) -> Optional[Tuple[str, str, str, str]]:
        """Extract barcode and primer info from golden string in sequence"""
        sequence = str(record.seq)

        # Look for golden string between delimiters
        start = sequence.find(START_DELIMITER)
        if start == -1:
            sequence = reverse_complement(sequence)
            start = sequence.find(START_DELIMITER)
            if start == -1:
                raise ValueError(format(f"Golden string does not include {START_DELIMITER} {sequence}"))
            else:
                logging.debug(f"Sequence was reversed after demultiplexing")

        end = sequence.find(END_DELIMITER, start)
        if end == -1:
            return None

        # Extract info between delimiters
        start_pos = start + len(START_DELIMITER)
        info_seq = sequence[start_pos:end]

        try:
            # Extract components using known lengths from Index.txt
            pos = 0

            # Extract forward barcode
            fwd_bc = info_seq[pos:pos + self.fwd_barcode_len]
            pos += self.fwd_barcode_len

            # Extract forward primer
            fwd_primer = info_seq[pos:pos + self.fwd_primer_len]
            pos += self.fwd_primer_len

            # Extract reverse primer
            rev_primer = info_seq[pos:pos + self.rev_primer_len]
            pos += self.rev_primer_len

            # Extract reverse barcode
            rev_bc = info_seq[pos:pos + self.rev_barcode_len]

            # Verify we got complete components
            if len(fwd_bc) != self.fwd_barcode_len or len(rev_bc) != self.rev_barcode_len:
                logging.error(f"Inconsistent barcode length for {info_seq} extracted from {sequence}")
                return None

            return (fwd_bc, fwd_primer, rev_primer, rev_bc)

        except (IndexError, ValueError):
            return None

    def validate_sample_file(self, sample_file: Path):
        """Validate reads in a single sample file"""
        # Extract specimen ID from filename
        sample_id = sample_file.stem.replace("sample_", "")

        is_fastq = sample_file.suffix.lower() in {'.fastq', '.fq'}
        format_type = 'fastq' if is_fastq else 'fasta'

        for record in SeqIO.parse(sample_file, format_type):
            self.global_stats.total_reads += 1

            # Always try to extract golden string info
            golden_info = self.extract_golden_info(record)
            if golden_info is None:
                self.global_stats.failed_golden_extraction += 1
                continue

            fwd_bc, fwd_primer, rev_primer, rev_bc = golden_info

            # Track barcode stats
            for bc in (fwd_bc, rev_bc):
                self.barcode_stats[bc].total_reads += 1

            # Handle unknown/ambiguous cases
            if sample_id.startswith(('unknown', 'ambiguous')):
                if sample_id.startswith('unknown'):
                    self.global_stats.unknown_reads += 1
                else:
                    self.global_stats.ambiguous_reads += 1

                # Track as unassigned for barcode stats
                for bc in (fwd_bc, rev_bc):
                    self.barcode_stats[bc].unassigned_reads += 1

                # Record edit distances if available
                if match := re.search(r'\((\d+),(\d+),(\d+),(\d+)\)', record.description):
                    _, b1_dist, _, b2_dist = map(int, match.groups())
                    self.barcode_stats[fwd_bc].edit_distances[b1_dist] += 1
                    self.barcode_stats[rev_bc].edit_distances[b2_dist] += 1
                continue

            # Check if specimen ID matches barcode pair
            correct_barcodes = False
            if sample_id in self.specimens:
                expected_spec = self.specimens[sample_id]
                if expected_spec.fwd_barcode == fwd_bc and expected_spec.rev_barcode == rev_bc:
                    correct_barcodes = True

            if correct_barcodes:
                self.global_stats.correct_assignments += 1
                for bc in (fwd_bc, rev_bc):
                    self.barcode_stats[bc].correct_reads += 1
            else:
                self.global_stats.incorrect_assignments += 1
                for bc in (fwd_bc, rev_bc):
                    self.barcode_stats[bc].incorrect_reads += 1

            # Record edit distances if available
            if match := re.search(r'\((\d+),(\d+),(\d+),(\d+)\)', record.description):
                _, b1_dist, _, b2_dist = map(int, match.groups())
                self.barcode_stats[fwd_bc].edit_distances[b1_dist] += 1
                self.barcode_stats[rev_bc].edit_distances[b2_dist] += 1

    def generate_report(self) -> str:
        """Generate validation report with fixed-width columns"""
        report = []

        # Overall statistics
        report.append("=== Overall Statistics ===")

        # Create fixed-width format for overall stats
        stat_format = "{:<25} {:>12,} {:>9.2f}%"

        total = self.global_stats.total_reads
        report.append(f"\nTotal Reads Processed: {total:,}")

        report.append("\nCategory Breakdown:")
        report.append("-" * 48)  # Header separator
        report.append(stat_format.format(
            "Correct",
            self.global_stats.correct_assignments,
            self.global_stats.accuracy()))
        report.append(stat_format.format(
            "Contaminant",
            self.global_stats.incorrect_assignments,
            self.global_stats.contamination_rate()))
        report.append(stat_format.format(
            "Unassigned",
            self.global_stats.unknown_reads + self.global_stats.ambiguous_reads,
            self.global_stats.unassigned_rate()))

        # Add failed extractions separately since they're not part of the main categories
        report.append(f"\nFailed Golden String Extraction: {self.global_stats.failed_golden_extraction:,}")

        # Barcode statistics
        report.append("\n\n=== Barcode Statistics ===")

        # Create fixed-width format for barcode stats
        header_format = "{:<12} {:>12} {:>12} {:>12} {:>12} {:>12}"
        data_format = "{:<12} {:>12,} {:>11.2f}% {:>11.2f}% {:>11.2f}% {:>11.2f}"

        report.append(header_format.format(
            "Barcode", "Total", "Correct", "Contaminant", "Unassigned", "Mean Edit"))
        report.append("-" * 75)  # Header separator

        # Sort barcodes by correct percentage in descending order
        sorted_barcodes = sorted(
            [(barcode, stats) for barcode, stats in self.barcode_stats.items() if stats.total_reads > 0],
            key=lambda x: x[1].accuracy(),
            reverse=True
        )

        for barcode, stats in sorted_barcodes:
            report.append(data_format.format(
                barcode,
                stats.total_reads,
                stats.accuracy(),
                stats.contamination_rate(),
                stats.unassigned_rate(),
                stats.mean_edit_distance()
            ))

        return "\n".join(report)

    def validate_demux_output(self, output_dir: Path):
        """Validate all sample files in output directory"""
        sample_files = list(output_dir.glob("sample_*.fa*"))
        sample_files.extend(output_dir.glob("sample_*.fq*"))

        if not sample_files:
            raise ValueError(f"No sample files found in {output_dir}")

        for sample_file in sample_files:
            self.validate_sample_file(sample_file)


def main():
    parser = argparse.ArgumentParser(
        description='Validate demultiplexing results using simulated reads')

    parser.add_argument('index_file', type=Path,
                        help='Index.txt file used for demultiplexing')
    parser.add_argument('demux_dir', type=Path,
                        help='Directory containing demultiplexed output files')
    parser.add_argument('--output', type=Path,
                        help='Output file for validation report (default: stdout)')

    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Initialize and run validator
    validator = DemuxValidator()
    logging.info("Reading index file...")
    validator.read_index_file(args.index_file)

    logging.info("Validating demultiplexed files...")
    validator.validate_demux_output(args.demux_dir)

    # Generate report
    report = validator.generate_report()

    if args.output:
        with open(args.output, 'w') as f:
            f.write(report)
        logging.info(f"Report written to {args.output}")
    else:
        print(report)


if __name__ == "__main__":
    main()

