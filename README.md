# Parallel Barcode Generator

A Python program that systematically generates DNA barcodes for dual-indexed sequencing, optimized for nanopore applications.

## Key Features

- **Exhaustive Search**: Unlike random sampling approaches, systematically evaluates all possible sequences to find optimal barcodes
- **Levenshtein Distance**: Uses edit distance instead of Hamming distance, providing better protection against insertions/deletions
- **Parallel Processing**: Utilizes multiple CPU cores for efficient generation of barcodes
- **Reverse Complement Checking**: Ensures barcodes are distinct from reverse complements of all other barcodes
- **Exclusion Sequences**: Can exclude sequences (primers, adapters, etc.) and their subsequences from consideration
- **GC Content Control**: Filters barcodes based on GC content range
- **Homopolymer Filtering**: Prevents long runs of identical bases

## Usage

```bash
parallel_barcode_generator.py [options]
```

### Options

- `--length`: Length of barcodes to generate (default: 13)
- `--number`: Number of barcodes to generate (default: 96)
- `--min-distance`: Minimum Levenshtein distance between barcodes (default: 7)
- `--min-gc`: Minimum GC content (0-1) (default: 0.4)
- `--max-gc`: Maximum GC content (0-1) (default: 0.6)
- `--max-homopolymer`: Maximum allowed homopolymer length (default: 4)
- `--exclusion-file`: FASTA file containing sequences to exclude
- `--processes`: Number of worker processes (default: CPU count)
- `--chunk-size`: Number of sequences to process per work item (default: 10000)
- `--output`: Output file path (default: barcode_<date>.txt)

## Output

The program produces a text file containing:
1. Generation parameters and statistics
2. List of generated barcodes with their GC content

## Attribution

This is a substantial revision of the barcode_generator program originally created by Luca Comai (Plant Biology and Genome Center, UC Davis). While it shares the same goal of generating DNA barcodes, it uses different algorithms and approaches.

## License

MIT License (same as original barcode_generator)