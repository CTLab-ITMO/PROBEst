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
                        help="Input FASTA file for probe generation. Probes are generated for different contigs separately. Only gene-coding regions are recommended (.fna).")

    parser.add_argument("-tb", "--true_base",
                        required=True,
                        help="Path to the BLAST database for primer adjustment. This database is used to ensure primer specificity.")

    parser.add_argument("-fb", "--false_base",
                        required=True,
                        nargs="*",
                        help="Path(s) to the BLAST database(s) for non-specific testing. These databases are used to filter out non-specific probes. Wildcards are not accepted.")

    parser.add_argument("-c", "--contig_table",
                        required=True,
                        help="Path to a .tsv table containing BLAST database information.")

    parser.add_argument("-o", "--output",
                        required=True,
                        help="Output directory path for storing results.")

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
                        help="Mutation probability per position in a primer. Default is 0.05.")

    parser.add_argument("-S", "--set_size",
                        required=False,
                        default=10, type=int,
                        help="Number of mutated probes to generate per primer. Default is 10.")

    parser.add_argument("-A", "--append",
                        required=False,
                        default=True, type=bool,
                        help="Whether to append the best probes to the array in the evolutionary algorithm. Default is True.")

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
    
    return args