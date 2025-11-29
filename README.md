# PROBEst v.0.1.4. <a href=""><img src="img/probest_logo.jpg" align="right" width="150" ></a> 
### St. Petersburg tool for genereting nucleotide probes with specified properties

[![python package](https://github.com/CTLab-ITMO/PROBEst/actions/workflows/build.yml/badge.svg)](https://github.com/CTLab-ITMO/PROBEst/actions/workflows/build.yaml)

<span style="color: red">**Warning**:</span> tool is under active development

**PROBEst** is a tool designed for generating nucleotide probes with specified properties, leveraging advanced algorithms and AI-driven techniques to ensure high-quality results. The tool is particularly useful for researchers and bioinformaticians who require probes with tailored universality and specificity for applications such as PCR, hybridization, and sequencing. By integrating a wrapped evolutionary algorithm, PROBEst optimizes probe generation through iterative refinement, ensuring that the final probes meet stringent biological and computational criteria.

At the core of PROBEst is an AI-enhanced workflow that combines Primer3 for initial oligonucleotide generation, BLASTn for specificity and universality checks, and a mutation module for probe optimization. The tool allows users to input target sequences, select reference files for universality and specificity validation, and customize layouts for probe design. The evolutionary algorithm iteratively refines the probes by introducing mutations and evaluating their performance, ensuring that the final output is both specific to the target and universally applicable across related sequences.


# Download and installation

## Installation

```bash
git clone https://github.com/CTLab-ITMO/PROBEst.git
cd PROBEst
conda env create -f environment.yml
conda activate probest
python setup.py install
```

### Validate installation

```bash
bash test_generator.sh
```

## Usage

### Preparation

`pipeline.py` relies on pre-prepared BLASTn databases. To create the required `true_base`, `false_base`, and `contig_table`, you can use the following script:

```bash
bash scripts/generator/prep_db.sh \
  -n {DATABASE_NAME} \
  -c {CONTIG_NAME} \
  -t {TMP_DIR} \
  [FASTA]
```

#### Arguments:
- `-n DATABASE_NAME`:  Name of the output BLAST database (required).  
- `-c CONTIG_TABLE`:  Output file to store contig names and their corresponding sequence headers (required).  
- `-t TMP_DIR`:  Temporary directory for intermediate files (optional, defaults to `./.tmp`).  
- `FASTA`:  List of input FASTA files (gzipped or uncompressed). 

### Generation

PROBEst can be run using the following command:

```bash
python pipeline.py \
  -i {INPUT} \
  -tb {TRUE_BASE} \
  -fb [FALSE_BASE ...] \
  -c {CONTIG_TABLE} \
  -o {OUTPUT}
```

**Blastn databases** and **contig table** are results of the ```prep_db.sh```

#### Key arguments:
- `-i INPUT`: Input FASTA file for probe generation.
- `-tb TRUE_BASE`: Input BLASTn database path for primer adjusting.
- `-fb FALSE_BASE`: Input BLASTn database path for non-specific testing.
- `-c CONTIG_TABLE`: .tsv table with BLAST database information.
- `-o OUTPUT`: Output path for results.
- `-t THREADS`: Number of threads to use.
- `-a ALGORITHM`: Algorithm for probe generation (`FISH` or `primer`).

For a full list of arguments, run:

```bash
python pipeline.py --help
```

For parameter selection, grid search is implemented. You can specify parameters in json (see for example `data/test/general/param_grid_light.json`) and run 

```bash
python test_parameters.py \
  -p {JSON}
```


# Algorithm

## Algorithm Steps

0. **Prepare BLASTn databases**

1. **Select File for Probe Generation** (`INPUT`)

2. **Select Files for Universality Check** (`TRUE_BASE`)

3. **Select Files for Specificity Check** (`FALSE_BASE`)
   
4. **Select Layouts and Run Wrapped Evolutionary Algorithm** (`pipeline.py`)
   
   a. **Primer3 Generation**
      
   b. **BLASTn Check**
      
   c. **Parsing**
      
   d. **Mutation in Probe**
   
   e. **AI corrections**

```mermaid
---
config:
  layout: elk
  look: classic
---
%%{init: {
  'theme': 'base',
  'themeVariables': {
    'fontFamily': 'arial',
    'fontSize': '16px',
    'primaryColor': '#fff',
    'primaryBorderColor': '#FFAC1C',
    'primaryTextColor': '#000',
    'lineColor': '#000',
    'secondaryColor': 'white',
    'tertiaryColor': '#fff',
    'subgraphBorderStyle': 'dotted'
  },
  'flowchart': {
    'curve': 'monotoneY',
    'padding': 15
  }
}}%%

graph LR
  subgraph inputs
  A
  A1
  T1
  T3
  end

  A([Initial probe generation]):::input -- primer3 --> B2(initial probe set):::probe
  A -- oligominer --> B2
  A1([Custom probes]):::input --> B2

  T1([Target sequences]):::input -- blastn-db --> T2[(target)]
  T3([Offtarget sequences]):::input -- blastn-db --> T4[(offtarget)]

  subgraph database
  T2
  T4
  end

  T2 --> EA
  T4 --> EA
  B2 --> EA

  EA[evolutionary algorithm] --> T11(results):::probe

  classDef empty width:0px,height:0px;
  classDef input fill:#90EE9020,stroke:#fff,stroke-width:2px,shape:ellipse;
  classDef probe fill:#FFAC1C20,stroke:#fff,stroke-width:2px;
```

```mermaid
---
config:
  layout: elk
  look: classic
---
%%{init: {
  'layout': 'elk',
  'theme': 'base',
  'themeVariables': {
    'fontFamily': 'arial',
    'fontSize': '16px',
    'primaryColor': '#fff',
    'primaryBorderColor': '#FFAC1C',
    'primaryTextColor': '#000',
    'lineColor': '#000',
    'secondaryColor': 'white',
    'tertiaryColor': '#fff',
    'subgraphBorderStyle': 'dotted'
  },
  'flowchart': {
    'curve': 'monotoneY',
    'padding': 15
  }
}}%%

graph LR
  subgraph evolutionary algorithm
    subgraph hits
      TP
      TN
    end

    B(probe set):::probe --> TP[target]
    B --> TN[offtarget]
    B1 -- mutations --> B

    TP -- coverage --> T6[universality]
    TP -- duplications --> T7[multimapping]
    TN ---> T8[specificity]

    subgraph check
    T6
    T7
    T8
    M1
    end

    B --- E6[ ]:::empty --> M1[modeling]
    TP --- E6

    M1 --- E3[ ]:::empty
    T6 --- E3
    T7 --- E3
    T8 --- E3
    E3 -- quality prediction --> B1(filtered probe set):::probe
  end
  B1 --> T11(results):::probe

  classDef empty width:0px,height:0px;
  classDef input fill:#90EE9020,stroke:#fff,stroke-width:2px,shape:ellipse;
  classDef probe fill:#FFAC1C20,stroke:#fff,stroke-width:2px;
```
    

## Project Structure

```mermaid
---
config:
  theme: neutral
  look: classic
---
%%{init: {
  'theme': 'base',
  'themeVariables': {
    'fontFamily': 'arial',
    'fontSize': '16px',
    'primaryColor': '#fff',
    'primaryBorderColor': '#FFAC1C',
    'primaryTextColor': '#000',
    'lineColor': '#000',
    'secondaryColor': '#90EE90',
    'tertiaryColor': '#fff',
    'subgraphBorderStyle': 'dotted'
  },
  'flowchart': {
    'curve': 'monotoneY',
    'padding': 15
  }
}}%%

graph LR
    PROBEst([PROBEst]) --> src[src/]
    PROBEst --> scripts[scripts/]
    PROBEst --> tests[tests/]

    subgraph folders
    src
    scripts
    tests
    end
    
    src --> C[benchmarking]
    src --> A[generation]
    tests --> A
    
    scripts --> D[preprocessing]
    scripts --> B[database parsing]
    D --> A
```

# Testing

- To check the installation: `bash test_generator.sh`

- For developers: use `pytest`


# License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

# Contribution

We welcome contributions from the community! To contribute:


Please read the [Contribution Guidelines](CONTRIBUTING.md) for more details.

# Wiki

Tool have its own <a href = "https://github.com/CTLab-ITMO/PROBEst/wiki">Wiki</a> pages with detailed information on usage cases, data description and another neccessary information
