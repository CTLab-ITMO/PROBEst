# MIT License
#
# Copyright (c) 2025 CTLab-ITMO
#
# Authors: Daniil Smutin, Aleksandr Serdiukov, Vitalii Dravgelis, Artem Ivanov,
# Aleksei Zabashta, Sergey Muravyov, and the CTLab-ITMO university team.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import argparse

def arguments_parse():
    """
    Parses command-line arguments for the probe generation tool.

    This function defines and parses the command-line arguments required for generating probes
    based on FASTA files and BLAST databases. It includes options for input files, output paths,
    algorithm settings, and parameters for primer3 and blastn.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.
    """
    
    
    description = """
    Generation of probes based on FASTA files and BLAST databases.

    PROBESt generates probes for gene-coding regions using an evolutionary algorithm.
    It requires primer3 and blastn to be pre-installed. The tool allows for the generation
    of initial primer sets, adjustment of primers based on BLAST databases, and removal of
    non-specific probes.

    Usage:
    1. Provide an input FASTA file for probe generation.
    2. Specify BLAST databases for primer adjustment and non-specific testing.
    3. Configure evolutionary algorithm parameters and primer3/blastn settings.
    """

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    # Main arguments
    parser.add_argument("-i", "--input",
                        required=True,
                        nargs="*",
                        help="Input FASTA file(s) or directory(ies) for probe generation. Can be a single file/directory or multiple files/directories. If a directory is provided, all *.fa, *.fna, *.fasta files (and their .gz versions) will be processed. Probes are generated for different contigs separately. Only gene-coding regions are recommended (.fna).")

    parser.add_argument("-tb", "--true_base",
                        required=True,
                        help="Path to the BLAST database for primer adjustment. This database is used to ensure primer specificity.")

    parser.add_argument("-fb", "--false_base",
                        required=True,
                        nargs="*",
                        help="Path(s) to the BLAST database(s) for non-specific testing. These databases are used to filter out non-specific probes. Wildcards are not accepted.")

    parser.add_argument("-c", "--contig_table",
                        required=False,
                        default=None,
                        help="Path to a .tsv table containing BLAST database information. If not provided and FASTA directories are used for bases, it will be auto-generated in the output directory.")

    parser.add_argument("-o", "--output",
                        required=True,
                        help="Output directory path for storing results.")

    # BLAST database preparation arguments (for FASTA directory inputs)
    parser.add_argument("--prep_db_tmp",
                        required=False,
                        default=None,
                        help="Temporary directory for BLAST database preparation when using FASTA directories as input. If not provided, a temporary directory will be created in the output directory.")

    parser.add_argument("-t", "--threads",
                        required=False,
                        default="1",
                        help="Number of threads to use for parallel processing. Default is 1.")

    parser.add_argument("-a", "--algorithm",
                        required=False,
                        default="FISH",
                        help="Algorithm for probe generation. Options: 'FISH' (default) or 'primer'.")

    parser.add_argument("-ot", "--output_tmp",
                        required=False,
                        default="",
                        help="Path to the temporary output directory for calculations and data processing. Default is a '.tmp' directory within the output directory.")

    # Evolutionary algorithm arguments
    parser.add_argument("-N", "--iterations",
                        required=False,
                        default=5, type=int,
                        help="Maximum number of iterations for the evolutionary algorithm. Default is 5.")

    parser.add_argument("-T", "--top",
                        required=False,
                        default=10, type=int,
                        help="Number of top probes to mutate and use in the next generation. Default is 10.")

    parser.add_argument("-M", "--mutation_rate",
                        required=False,
                        default=0.05, type=float,
                        help="'SNP' mutation probability per position in a primer. Default is 0.05.")
    
    parser.add_argument("-I", "--indel_rate",
                        required=False,
                        default=0.05, type=float,
                        help="'InDel' mutation probability per position in a primer. Default is 0.05.")

    parser.add_argument("-S", "--set_size",
                        required=False,
                        default=10, type=int,
                        help="Number of mutated probes to generate per primer. Default is 10.")

    parser.add_argument("-A", "--append",
                        required=False,
                        default=True, type=bool,
                        help="Whether to append the best probes to the array in the evolutionary algorithm. Default is True.")

    # De-degeneration algorithm arguments
    parser.add_argument("-D", "--dedegeneration_iterations",
                        required=False,
                        default=5, type=int,
                        help="Maximum number of iterations for the de-degeneration algorithm. Default is 5.")
    
    parser.add_argument("--dedegeneration_set_size",
                        required=False,
                        default=10, type=int,
                        help="Number of de-degenerated variants to generate per probe in the de-degeneration algorithm. Default is 10.")
    
    parser.add_argument("--dedegeneration_append",
                        required=False,
                        default=True, type=bool,
                        help="Whether to append the best probes to the array in the de-degeneration algorithm. Default is True.")
    
    parser.add_argument("--dedegeneration_rate",
                        required=False,
                        default=0.5, type=float,
                        help="De-degeneration probability per degenerate nucleotide position. Default is 0.5.")

    # Execution arguments
    parser.add_argument("--primer3",
                        required=False,
                        default="primer3_core",
                        help="Path or command to execute primer3_core. Default is 'primer3_core'.")

    parser.add_argument("--blastn",
                        required=False,
                        default="blastn",
                        help="Path or command to execute blastn. Default is 'blastn'.")

    parser.add_argument("--add_set",
                        required=False,
                        default=None,
                        nargs="*",
                        help="File(s) containing additional sets of probes to append to the initial primer3 generation. Default is None.")

    # Initial set generator selection
    parser.add_argument("--initial_generator",
                        required=False,
                        default="primer3",
                        choices=["primer3", "oligominer"],
                        help="Tool to use for initial probe set generation. Options: 'primer3' (default) or 'oligominer'.")

    # OligoMiner arguments
    parser.add_argument("--oligominer_path",
                        required=False,
                        default=None,
                        help="Path to the OligoMiner installation directory. Required when --initial_generator=oligominer.")

    parser.add_argument("--oligominer_python",
                        required=False,
                        default=None,
                        help="Python interpreter to use for OligoMiner scripts (e.g., 'python2', 'python2.7', or path to Python 2). If not specified, will try to auto-detect Python 2, then fall back to 'python'.")

    parser.add_argument("--oligominer_probe_length",
                        required=False,
                        default=None,
                        type=int,
                        help="OligoMiner option: Probe length in nucleotides. If not specified, blockParse.py will use its default (36-41).")

    parser.add_argument("--oligominer_temperature",
                        required=False,
                        default=None,
                        type=int,
                        help="OligoMiner option: Melting temperature in Celsius. If not specified, blockParse.py will use its default (42-47).")

    parser.add_argument("--oligominer_insert_coords",
                        required=False,
                        default=None,
                        help="OligoMiner option: Path to BED file with insert coordinates for filtering probes. If not provided, all probes are kept.")

    parser.add_argument("--oligominer_keep_tmp",
                        required=False,
                        default=False,
                        type=bool,
                        help="OligoMiner option: Whether to keep temporary files after processing. Default is False.")

    # Primer3 template arguments
    parser.add_argument("--PRIMER_PICK_PRIMER",
                        required=False,
                        default=10,
                        help="Primer3 template option: Number of probes to pick. Default is 10.")

    parser.add_argument("--PRIMER_NUM_RETURN",
                        required=False,
                        default=10,
                        help="Primer3 template option: Initial set size per gene. Default is 10.")

    parser.add_argument("--PRIMER_OPT_SIZE",
                        required=False,
                        default=25,
                        type=int,
                        help="Primer3 template option: Optimal primer length. Default is 25.")

    parser.add_argument("--PRIMER_MIN_SIZE",
                        required=False,
                        default=15,
                        type=int,
                        help="Primer3 template option: Minimum primer length. Default is 15.")

    parser.add_argument("--PRIMER_MAX_SIZE",
                        required=False,
                        default=30,
                        type=int,
                        help="Primer3 template option: Maximum primer length. Default is 30.")

    parser.add_argument("--PRIMER_PRODUCT_SIZE_RANGE",
                        required=False,
                        default="100-1000",
                        help="Primer3 template option: Desired PCR product size range. Provide two values separated by '-'. Default is '100-1000'.")

    # BLAST template arguments
    parser.add_argument("--blastn_base_prebuilt",
                        required=False,
                        help="Specify if the BLAST database is prebuilt. Default is None.")

    parser.add_argument("--word_size",
                        required=False,
                        default="7",
                        help="BLAST template option: Word size for alignment. Default is 7.")

    parser.add_argument("--reward",
                        required=False,
                        default="3",
                        help="BLAST template option: Reward for a match. Default is 3.")

    parser.add_argument("--penalty",
                        required=False,
                        default="-3",
                        help="BLAST template option: Penalty for a mismatch. Default is -3.")

    parser.add_argument("--gapopen",
                        required=False,
                        default="6",
                        help="BLAST template option: Penalty for opening a gap. Default is 6.")

    parser.add_argument("--gapextend",
                        required=False,
                        default="3",
                        help="BLAST template option: Penalty for extending a gap. Default is 3.")

    parser.add_argument("--evalue",
                        required=False,
                        default="1",
                        help="BLAST template option: E-value threshold. Default is 1.")

    # Probe check template arguments
    parser.add_argument("--max_mismatch",
                        required=False,
                        default="5",
                        help="Probe check template option: Maximum allowed mismatches. Default is 5.")

    parser.add_argument("--multimap_max",
                        required=False,
                        default="1",
                        help="Probe check template option: Maximum allowed multimapped hits. Default is 1.")

    parser.add_argument("--negative_max",
                        required=False,
                        default="0",
                        help="Probe check template option: Maximum allowed negative hits. Default is 0.")

    parser.add_argument("--min_ident",
                        required=False,
                        default="70",
                        help="Probe check template option: Minimum identity percentage. Default is 70.")

    parser.add_argument("--visualize",
                        required=False,
                        default=False,
                        type=bool,
                        help="Whether to create visualizations for probe-target pairs. Default is False.")

    parser.add_argument("--AI",
                        required=False,
                        default=True,
                        type=bool,
                        help="Whether to apply AI-based filtration to probes. Default is True. Set to False to disable AI filtration.")

    args = parser.parse_args()
    
    # Correct args
    # Correct true/false bases. Trim trailing spaces in db path to avoid blastn crash
    args.true_base = args.true_base[:-
                                    1] if args.true_base.endswith('/') else args.true_base
    for i, db_neg in enumerate(args.false_base):
        if db_neg.endswith('/'):
            args.false_base[i] = db_neg[:-1]
            
    # Correct tmp
    if args.output_tmp == "":
        args.output_tmp = args.output + "/.tmp/"
    
    # Set default contig_table if not provided
    if args.contig_table is None:
        args.contig_table = args.output + "/contigs.tsv"
    
    return args