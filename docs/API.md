# PROBEst API Documentation

## Overview

PROBEst (Probe Generation and Optimization Tool) is a Python package designed for generating and optimizing nucleotide probes with specified properties. The tool combines Primer3 for initial oligonucleotide generation, BLASTn for specificity and universality checks, and an evolutionary algorithm for probe optimization.

## Core Components

### 1. Command Line Interface (`args.py`)

The main entry point for the PROBEst tool is through its command-line interface, which provides extensive configuration options for probe generation and optimization.

#### Main Arguments

- `-i, --input`: Input FASTA file for probe generation (required)
- `-tb, --true_base`: BLAST database path for primer adjustment (required)
- `-fb, --false_base`: BLAST database path(s) for non-specific testing (required)
- `-c, --contig_table`: Path to .tsv table with BLAST database information (required)
- `-o, --output`: Output directory path (required)
- `-t, --threads`: Number of threads for parallel processing (default: 1)
- `-a, --algorithm`: Probe generation algorithm ('FISH' or 'primer', default: 'FISH')

#### Evolutionary Algorithm Parameters

- `-N, --iterations`: Maximum iterations (default: 5)
- `-T, --top`: Number of top probes to mutate (default: 10)
- `-M, --mutation_rate`: SNP mutation probability (default: 0.05)
- `-I, --indel_rate`: InDel mutation probability (default: 0.05)
- `-S, --set_size`: Mutated probes per primer (default: 10)
- `-A, --append`: Append best probes to array (default: True)

#### Primer3 Configuration

- `--PRIMER_PICK_PRIMER`: Number of probes to pick (default: 10)
- `--PRIMER_NUM_RETURN`: Initial set size per gene (default: 10)
- `--PRIMER_OPT_SIZE`: Optimal primer length (default: 25)
- `--PRIMER_MIN_SIZE`: Minimum primer length (default: 15)
- `--PRIMER_MAX_SIZE`: Maximum primer length (default: 30)
- `--PRIMER_PRODUCT_SIZE_RANGE`: PCR product size range (default: "100-1000")

#### BLAST Configuration

- `--word_size`: Word size for alignment (default: "7")
- `--reward`: Match reward (default: "3")
- `--penalty`: Mismatch penalty (default: "-3")
- `--gapopen`: Gap opening penalty (default: "6")
- `--gapextend`: Gap extension penalty (default: "3")
- `--evalue`: E-value threshold (default: "1")

### 2. Evolution Module (`evolution.py`)

The evolution module implements the core mutation and optimization algorithms for probe refinement.

#### Functions

##### `mutate_position(x, mutrate, indelrate)`

Mutates a nucleotide based on given mutation and indel rates.

**Parameters:**
- `x` (str): Nucleotide to potentially mutate
- `mutrate` (float): Mutation rate
- `indelrate` (float): Rate of insertions/deletions

**Returns:**
- `str`: Mutated nucleotide or empty string (for deletions)

##### `mutate_sequence(args, out_dir, seqs)`

Applies mutations to a set of sequences and writes results to a FASTA file.

**Parameters:**
- `args`: Command line arguments
- `out_dir` (str): Output directory path
- `seqs` (dict): Dictionary of sequences to mutate

### 3. Genome Operations (`genome_operations.py`)

Handles genome-related operations including fetching sequences and performing BLAST searches.

### 4. Probe Alignment Profiler (`probe_alignment_profiler.py`)

Analyzes and profiles probe alignments for specificity and universality.

### 5. Primer3 Integration (`primer3.py`)

Manages the integration with Primer3 for initial probe generation.

### 6. Merge Module (`merge.py`)

Handles merging and combining probe sets from different sources.

### 7. Miscellaneous Utilities (`misc.py`)

Provides utility functions for:
- Writing statistics
- FASTA file operations
- Output directory management
- Sequence pairing

## Usage Examples

### Basic Probe Generation

```bash
python pipeline.py \
    -i input.fasta \
    -tb true_base_db \
    -fb false_base_db1 false_base_db2 \
    -c contig_table.tsv \
    -o output_directory
```

### Advanced Configuration

```bash
python pipeline.py \
    -i input.fasta \
    -tb true_base_db \
    -fb false_base_db1 false_base_db2 \
    -c contig_table.tsv \
    -o output_directory \
    -t 4 \
    -a FISH \
    -N 10 \
    -T 20 \
    -M 0.1 \
    -I 0.05 \
    -S 15
```

## Dependencies

- Primer3
- BLASTn
- Python 3.6+
- Required Python packages (see requirements.txt)

## Notes

- The tool is under active development
- Some features may require specific versions of dependencies
- Performance may vary based on input data size and system resources 