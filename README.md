# PROBEst v.0.1. <a href=""><img src="img/probest_logo.jpg" align="right" width="150" ></a> 
### St. Petersburg tool for genereting nucleotide probes with specified properties

PROBEst is a sophisticated tool designed for generating nucleotide probes with specified properties, leveraging advanced algorithms and AI-driven techniques to ensure high-quality results. The tool is particularly useful for researchers and bioinformaticians who require probes with tailored universality and specificity for applications such as PCR, hybridization, and sequencing. By integrating a wrapped evolutionary algorithm, PROBEst optimizes probe generation through iterative refinement, ensuring that the final probes meet stringent biological and computational criteria.

At the core of PROBEst is an AI-enhanced workflow that combines Primer3 for initial primer generation, BLASTn for specificity and universality checks, and a mutation module for probe optimization. The tool allows users to input target sequences, select reference files for universality and specificity validation, and customize layouts for probe design. The evolutionary algorithm iteratively refines the probes by introducing mutations and evaluating their performance, ensuring that the final output is both specific to the target and universally applicable across related sequences. This AI-driven approach significantly enhances the efficiency and accuracy of probe generation, making PROBEst a valuable resource for molecular biology research.

**Warning**: tool is under development

# Algorithm

This section outlines the algorithm used for probe generation, including checks for universality and specificity.

## Algorithm Steps

1. **Select File for Probe Generation**
   - Choose the primary file that will be used to generate the probe.

2. **Select Files for Universality Check**
   - Identify and select files that will be used to assess the universality of the probe.

3. **Select Files for Specificity Check**
   - Identify and select files that will be used to evaluate the specificity of the probe.

4. **Select Layouts**
   - Determine the layouts that will be utilized during the probe generation process.

5. **Run Wrapped Evolutionary Algorithm**
   - Execute the following steps within the evolutionary algorithm:
   
   a. **Primer3 Generation**
      - Generate primers using the Primer3 tool.
      
   b. **BLASTn Check**
      - Perform a BLASTn check to ensure the generated probes are suitable.
      
   c. **Parsing**
      - Parse the results from the BLASTn check to extract relevant information.
      
   d. **Mutation in Probe**
      - Introduce mutations in the probe based on the parsed data to optimize performance.
  
   e. **AI corrections**
      - Probe evaluation based on AI optimizing function


# Download and installation

```bash
git clone github.com/CTLab-ITMO/PROBEst/
python setup.py develop
```

# Testing

- To check the installation: `bash test_run_generator.sh`

- For developers: use `.test/`

# Main usage

```
usage: pipeline.py [-h] -i INPUT -tb TRUE_BASE -fb [FALSE_BASE [FALSE_BASE ...]] -c
                   CONTIG_TABLE -o OUTPUT [-t THREADS] [-a ALGORITHM] [-ot OUTPUT_TMP]
                   [-N ITERATIONS] [-T TOP] [-M MUTATION_RATE] [-S SET_SIZE] [-A APPEND]
                   [--primer3 PRIMER3] [--blastn BLASTN] [--add_set [ADD_SET [ADD_SET ...]]]
                   [--PRIMER_PICK_PRIMER PRIMER_PICK_PRIMER]
                   [--PRIMER_NUM_RETURN PRIMER_NUM_RETURN]
                   [--PRIMER_OPT_SIZE PRIMER_OPT_SIZE] [--PRIMER_MIN_SIZE PRIMER_MIN_SIZE]
                   [--PRIMER_MAX_SIZE PRIMER_MAX_SIZE]
                   [--PRIMER_PRODUCT_SIZE_RANGE PRIMER_PRODUCT_SIZE_RANGE]
                   [--word_size WORD_SIZE] [--reward REWARD] [--penalty PENALTY]
                   [--gapopen GAPOPEN] [--gapextend GAPEXTEND] [--evalue EVALUE]
                   [--max_mismatch MAX_MISMATCH] [--multimap_max MULTIMAP_MAX]
                   [--negative_max NEGATIVE_MAX] [--min_ident MIN_IDENT]

Generation of probes based on fasta-files and blastn databases. To use it, select one
reference file to generate the initial primer set; blastn base to check primer universality
and cut off multimapping; blastn bases to remove non-specific probes Requires primer3 and
blastn pre-installed

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input FASTA file for generation. probes are generated for different
                        contigs separatly. Only gene-coding regions recommended (.fna)
  -tb TRUE_BASE, --true_base TRUE_BASE
                        Input blastn database path for primer adjusting
  -fb [FALSE_BASE [FALSE_BASE ...]], --false_base [FALSE_BASE [FALSE_BASE ...]]
                        Input blastn database path for non-specific testing. Wildcards are
                        not accepted
  -c CONTIG_TABLE, --contig_table CONTIG_TABLE
                        .tsv table with blast db information
  -o OUTPUT, --output OUTPUT
                        Output path
  -t THREADS, --threads THREADS
                        number of threads
  -a ALGORITHM, --algorithm ALGORITHM
                        algorithm for probes generation. 'FISH' as default, also could be
                        'primer'
  -ot OUTPUT_TMP, --output_tmp OUTPUT_TMP
                        Output .tmp dicrectory path for calculations and data processing.
                        .tmp in output directory as default
  -N ITERATIONS, --iterations ITERATIONS
                        Maximum iterations of evolutionary algorithm. 100 by default
  -T TOP, --top TOP     Top probes to mutate and use in next generation
  -M MUTATION_RATE, --mutation_rate MUTATION_RATE
                        Mutation probability per position of primer
  -S SET_SIZE, --set_size SET_SIZE
                        Size of mutated probes per primer
  -A APPEND, --append APPEND
                        Append best probes to array in evolutionary algoritm
  --primer3 PRIMER3     primer3_core path or command to exec. 'primer3' as default
  --blastn BLASTN       blastn path or command to exec. 'blastn' as default
  --add_set [ADD_SET [ADD_SET ...]]
                        file to set of probes to append to initial primer3 generation. empty
                        by default
  --PRIMER_PICK_PRIMER PRIMER_PICK_PRIMER
                        primer3 template option. Number of probes to pick
  --PRIMER_NUM_RETURN PRIMER_NUM_RETURN
                        primer3 template option. initial set size per gene
  --PRIMER_OPT_SIZE PRIMER_OPT_SIZE
                        primer3 template option
  --PRIMER_MIN_SIZE PRIMER_MIN_SIZE
                        primer3 template option
  --PRIMER_MAX_SIZE PRIMER_MAX_SIZE
                        primer3 template option
  --PRIMER_PRODUCT_SIZE_RANGE PRIMER_PRODUCT_SIZE_RANGE
                        primer3 template option. 2 values sepatated by '-'
  --word_size WORD_SIZE
                        blastn template option
  --reward REWARD       blastn template option
  --penalty PENALTY     blastn template option
  --gapopen GAPOPEN     blastn template option
  --gapextend GAPEXTEND
                        blastn template option
  --evalue EVALUE       blastn template option
  --max_mismatch MAX_MISMATCH
                        probe_check template option. maximum avialable mismatch
  --multimap_max MULTIMAP_MAX
                        probe_check template option. maximum multimapped hits
  --negative_max NEGATIVE_MAX
                        probe_check template option. maximum negative hits
  --min_ident MIN_IDENT
                        probe_check template option. minimal identity, percent
```
```
```
