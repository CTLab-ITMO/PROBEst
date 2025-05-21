# PROBEst generation

## Overview

The PROBEst probe generation tool (`pipeline.py`) is the main executable script for generating and optimizing nucleotide probes. It implements a sophisticated workflow that combines initial probe generation, evolutionary optimization, and quality control steps.

## Workflow

The pipeline follows these main steps:

1. **Initial Set Generation**
   - Parses input FASTA file
   - Generates initial probe set using Primer3
   - Merges and processes initial probes

2. **Evolutionary Algorithm Loop**
   - Performs BLASTn searches against reference databases
   - Detects and filters multimapping probes
   - Matches and evaluates probe performance
   - Applies mutations for optimization (if not final iteration)

3. **Output Generation**
   - Creates final FASTA file with optimized probes
   - Generates statistics report

## Usage

### Basic Usage

```bash
python pipeline.py \
    -i input.fasta \
    -tb true_base_db \
    -fb false_base_db1 false_base_db2 \
    -c contig_table.tsv \
    -o output_directory
```

### Advanced Usage

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

### Parameter Grid Search

PROBEst includes a grid search tool (`test_parameters.py`) for evaluating different parameter combinations. This is useful for:
- Finding optimal parameter settings
- Evaluating parameter sensitivity
- Comparing different algorithm configurations

#### Grid Search Usage

1. **Create Parameter Grid File**
   Create a JSON file defining parameter ranges to test. Example:
   ```json
   {
       "mutation_rate": [0.05, 0.1, 0.15],
       "indel_rate": [0.05, 0.1],
       "set_size": [10, 15, 20],
       "iterations": [5, 10]
   }
   ```

2. **Run Grid Search**
   ```bash
   python test_parameters.py \
       -p param_grid.json \
       -o ./results/grid_search \
       -t 8 \
       -tpw 2 \
       -bs 3
   ```

#### Grid Search Options

- `-p, --params_grid`: JSON file with parameter grid (required)
- `-o, --output`: Output directory (default: ./data/param_search)
- `-t, --threads`: Total number of threads
- `-tpw, --threads_per_worker`: Threads per pipeline instance
- `-kf, --keep_failed`: Include failed runs in results
- `-bs, --bootstrap`: Number of repetitions per combination

#### Output Analysis

The grid search generates:
- `param_stats.csv`: Comprehensive results table
- Individual run directories with:
  - Pipeline output files
  - Performance statistics
  - Log files

#### Best Practices

1. **Parameter Selection**
   - Start with wide parameter ranges
   - Focus on most influential parameters
   - Consider computational resources

2. **Resource Management**
   - Adjust thread count based on system
   - Monitor disk space for results
   - Use appropriate bootstrap count

3. **Analysis**
   - Compare max and mean hits
   - Evaluate convergence patterns
   - Consider parameter interactions

## Pipeline Components

### 1. Initial Set Generation

The pipeline starts by:
- Creating a temporary directory structure
- Converting input FASTA to single-line format
- Generating initial probe set using Primer3
- Merging and processing initial probes

### 2. Evolutionary Algorithm

For each iteration:
1. **BLASTn Analysis**
   - Searches against true base database
   - Searches against false base databases
   - Generates hit statistics

2. **Probe Filtering**
   - Detects multimapping probes
   - Filters based on specificity criteria
   - Calculates hit statistics

3. **Probe Matching**
   - Matches probes with their performance metrics
   - Calculates maximum and mean hits
   - Updates statistics

4. **Mutation (if not final iteration)**
   - Applies mutations to selected probes
   - Generates new probe variants
   - Merges and processes mutated probes

### 3. Output Generation

The pipeline produces:
- Final FASTA file with optimized probes
- Statistics report with performance metrics
- Temporary files for debugging (if needed)

## Output Files

### Main Output
- `output.fa`: Final FASTA file containing optimized probes
- `stats.json`: Statistics report with performance metrics

### Temporary Files (in output_tmp directory)
- `{iteration}/output.fa`: Probes for each iteration
- `{iteration}/merged.fa`: Merged probes
- `{iteration}/positive_hits.tsv`: BLAST hits against true base
- `{iteration}/negative_hits.tsv`: BLAST hits against false bases
- `{iteration}/clear_hits.tsv`: Filtered probe hits
- `{iteration}/probe_check/`: Probe check results

## Error Handling

The pipeline includes error handling for:
- Empty files after filtration
- Invalid probe configurations
- BLAST search failures
- File system issues

## Performance Considerations

- Use `-t` parameter to utilize multiple threads
- Adjust `-N` and `-S` parameters based on available memory
- Consider using `-A false` for memory-intensive operations
- Monitor temporary directory size during long runs

## Tips and Best Practices

1. **Input Preparation**
   - Ensure input FASTA is properly formatted
   - Verify BLAST databases are up-to-date
   - Check contig table format

2. **Parameter Selection**
   - Start with default parameters
   - Adjust mutation rates based on results
   - Fine-tune BLAST parameters for your use case

3. **Resource Management**
   - Monitor disk space for temporary files
   - Adjust thread count based on system resources
   - Consider using a dedicated output directory

4. **Quality Control**
   - Review statistics after each run
   - Check probe specificity
   - Validate results against known sequences 