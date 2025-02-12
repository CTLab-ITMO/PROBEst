# Pipeline: PROBESt----

# 1. Initial set generation
#
# < evolutionary algorithm >
# 2. blastn
# 3. multimapping detection and filtering
# 4. matching (depr)
# 5*. mutations (if not final)
# </ evolutionary algorithm >
#
# 6. output

# 0. Imports: python packages ----

import os
import argparse
from Bio import SeqIO
import subprocess
import numpy as np
import pandas as pd
import re
import random
import sys


# Add 'src' folder to sys.path for imports to work

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'src')))

# 0 Imports: absolute import ----
from PROBESt.primer3 import primer_template, parse_primer3_output
from PROBESt.evolution import mutate_seq
from PROBESt.misc import write_fasta, pairing  # out_dir
from PROBESt.merge import merge


def out_dir(iter):
    if args.output_tmp == "":
        return args.output + "/.tmp/" + str(iter) + "/"
    else:
        return args.output_tmp + str(iter) + "/"

# Merge left and right


def merge_iter(iter):
    out = out_dir(iter)
    merge(algo=args.algorithm,
            input=out + "output.fa",
            output=out + "merged.fa",
            tmp=out + "fasta_table.tsv",
            N=10,
            script_path=script_path)


def write_stats(stats: dict, output_dir: str):
    stats_file = output_dir + '/stats.csv'
    with open(stats_file, 'w') as f:
        f.write('iteration,max_hits,mean_hits\n')
        for iter, stats in stats.items():
            f.write(f'{iter},{stats['max_hits']},{stats['mean_hits']}\n')


# 0. Argparsing ----
description = "Generation of probes based on fasta-files and blastn databases.\n\nTo use it, select one reference file to generate the initial primer set; blastn base to check primer universality and cut off multimapping; blastn bases to remove non-specific probes\n\nRequires primer3 and blastn pre-installed"

parser = argparse.ArgumentParser(description=description)

# Main
parser.add_argument("-i", "--input",
                    required=True,
                    help="Input FASTA file for generation. probes are generated for different contigs separatly. Only gene-coding regions recommended (.fna)")

parser.add_argument("-tb", "--true_base",
                    required=True,
                    help="Input blastn database path for primer adjusting")

parser.add_argument("-fb", "--false_base",
                    required=True,
                    nargs="*",
                    help="Input blastn database path for non-specific testing. Wildcards are not accepted")

parser.add_argument("-c", "--contig_table",
                    required=True,
                    help=".tsv table with blast db information")

parser.add_argument("-o", "--output",
                    required=True,
                    help="Output path")

parser.add_argument("-t", "--threads",
                    required=False,
                    default="1",
                    help="number of threads")

parser.add_argument("-a", "--algorithm",
                    required=False,
                    default="FISH",
                    help="algorithm for probes generation. 'FISH' as default, also could be 'primer'")

parser.add_argument("-ot", "--output_tmp",
                    default="",
                    help="Output .tmp dicrectory path for calculations and data processing. .tmp in output directory as default")

# Evolutionary algoritm
parser.add_argument("-N", "--iterations",
                    default=5, type=int,
                    help="Maximum iterations of evolutionary algorithm. 100 by default")

parser.add_argument("-T", "--top",
                    default=10, type=int,
                    help="Top probes to mutate and use in next generation")

parser.add_argument("-M", "--mutation_rate",
                    default=0.05, type=float,
                    help="Mutation probability per position of primer")

parser.add_argument("-S", "--set_size",
                    default=10, type=int,
                    help="Size of mutated probes per primer")

parser.add_argument("-A", "--append",
                    default=True, type=bool,
                    help="Append best probes to array in evolutionary algoritm")

# Exec
parser.add_argument("--primer3",
                    required=False,
                    default="primer3",
                    help="primer3_core path or command to exec. 'primer3' as default")

parser.add_argument("--blastn",
                    required=False,
                    default="blastn",
                    help="blastn path or command to exec. 'blastn' as default")

parser.add_argument("--add_set",
                    required=False,
                    default=None,
                    nargs="*",
                    help="file to set of probes to append to initial primer3 generation. empty by default")

# Primer3 template
parser.add_argument("--PRIMER_PICK_PRIMER",
                    default=10,
                    help="primer3 template option. Number of probes to pick")

parser.add_argument("--PRIMER_NUM_RETURN",
                    default=10,
                    help="primer3 template option. initial set size per gene")

parser.add_argument("--PRIMER_OPT_SIZE",
                    default=25,
                    type=int,
                    help="primer3 template option")

parser.add_argument("--PRIMER_MIN_SIZE",
                    default=15,
                    type=int,
                    help="primer3 template option")

parser.add_argument("--PRIMER_MAX_SIZE",
                    default=30,
                    type=int,
                    help="primer3 template option")

parser.add_argument("--PRIMER_PRODUCT_SIZE_RANGE",
                    default="100-1000",
                    help="primer3 template option. 2 values sepatated by '-'")

# Blastn template
parser.add_argument("--word_size",
                    default="7",
                    help="blastn template option")

parser.add_argument("--reward",
                    default="3",
                    help="blastn template option")

parser.add_argument("--penalty",
                    default="-3",
                    help="blastn template option")

parser.add_argument("--gapopen",
                    default="6",
                    help="blastn template option")

parser.add_argument("--gapextend",
                    default="3",
                    help="blastn template option")

parser.add_argument("--evalue",
                    default="1",
                    help="blastn template option")

# probe_check template
parser.add_argument("--max_mismatch",
                    default="5",
                    help="probe_check template option. maximum avialable mismatch")

parser.add_argument("--multimap_max",
                    default="1",
                    help="probe_check template option. maximum multimapped hits")

parser.add_argument("--negative_max",
                    default="0",
                    help="probe_check template option. maximum negative hits")

parser.add_argument("--min_ident",
                    default="70",
                    help="probe_check template option. minimal identity, percent")

args = parser.parse_args()


# 1. Initial set generation ----
print("\n---- PROBESt v.0.1 ----\n")
print("Arguments passed")

script_path = os.path.dirname(
    os.path.realpath(__file__)) + "/scripts/generator/"
os.makedirs(out_dir(0), exist_ok=True)

# Make uniline fasta
uniline = "bash " + script_path + "uniline_fa.sh"
uniline += " -i " + args.input
uniline += " -o " + out_dir(0) + "input.fa"

subprocess.run(uniline, shell=True)
print("Input fasta parsed")

# Template generation
primer_temp = primer_template(
    out_dir(0) + "input.fa",
    args.PRIMER_PICK_PRIMER,
    args.PRIMER_OPT_SIZE,
    args.PRIMER_MIN_SIZE,
    args.PRIMER_MAX_SIZE,
    args.PRIMER_PRODUCT_SIZE_RANGE,
    args.PRIMER_NUM_RETURN)

template = open(out_dir(0)+"template", "w")
template.writelines(primer_temp)
template.close()

# Primer3 exec
primer3 = args.primer3 + " " + \
    out_dir(0) + "template" + " --output " + out_dir(0) + "output.p3"

subprocess.run(primer3, shell=True, executable="/bin/bash")
print("Primer3 done")

# Parse to fasta
probes = parse_primer3_output(out_dir(0) + "output.p3")
fasta = open(out_dir(0) + "output.fa", "w")
for primer in probes:
    sequence_id, primer_num, side, sequence = primer
    header = f">{sequence_id}_{primer_num}_{side}"
    fasta.write(f"{header}\n{sequence}\n")
fasta.close()

# Add probes to fasta
if args.add_set is not None:
    add_fasta = "cat " + args.add_set + " >> " + out_dir(0) + "output.fa"
    subprocess.run(add_fasta, shell=True, executable="/bin/bash")


merge_iter(0)

# < evolutionary algorithm >

# blastn command
blastn = args.blastn + " -num_threads " + \
    args.threads + " -outfmt '6 qseqid sseqid evalue sstart send ppos mismatch' " + \
    " -word_size " + args.word_size + \
    " -reward " + args.reward + \
    " -penalty " + args.penalty + \
    " -gapopen " + args.gapopen + \
    " -gapextend " + args.gapextend + \
    " -evalue " + args.evalue

# probe_check command
probe_check = "bash " + script_path + "/probe_check.sh" + \
    " -p " + script_path + "/probe_filt.py" + \
    " -d " + args.contig_table + \
    " -m " + str(args.top) + \
    " -e " + args.max_mismatch + \
    " -i " + args.min_ident + \
    " -a " + args.multimap_max + \
    " -b " + args.negative_max

stats = {}  # Hit stats for each iteration

for iter in range(1, args.iterations+1):
    print("\nIteration", iter, "----")
    os.makedirs(out_dir(iter), exist_ok=True)

    # 2. blastn ----
    blastn_iter = blastn + " -query " + out_dir(iter-1) + "merged.fa"

    # true base
    blastn_db = blastn_iter + " -db " + args.true_base + \
        " > " + out_dir(iter) + "positive_hits.tsv"
    subprocess.run(blastn_db, shell=True)

    print("Positive hits counted")

    # false bases
    for db_neg in args.false_base:
        blastn_db = blastn_iter + " -db " + db_neg + \
            " >> " + out_dir(iter) + "negative_hits.tsv"
        subprocess.run(blastn_db, shell=True)

    print("Negative hits counted")

    # 3. multimapping detection and filtering ----
    probe_check_iter = probe_check + \
        " -o " + out_dir(iter) + "clear_hits.tsv" + \
        " -r " + out_dir(iter) + "probe_check/" + \
        " -t " + out_dir(iter) + "positive_hits.tsv " + \
        out_dir(iter) + "negative_hits.tsv"

    subprocess.run(probe_check_iter, shell=True)

    # 4. probes matching ----
    try:
        probe_out = pd.read_table(out_dir(iter) + "clear_hits.tsv",
                                sep=' ', header=None)
    except:
        raise InterruptedError(
            "Empty file after filtration, try to use other probe_check properties and review false databases")

    probe_vals = probe_out.iloc[:, 0].value_counts()

    probe_list = list(set(probe_out.iloc[:, 0]))
    probe_list_hash = [hash(_) for _ in probe_list]

    max_hits = probe_vals.iloc[0]
    mean_hits = round(sum(probe_vals)/len(probe_vals), 1)
    stats[iter] = {
        'max_hits': max_hits,
        'mean_hits': mean_hits
    }
    print("Maximum hits:", max_hits)
    print("Mean hits:", mean_hits)

    # grep in probes.fa from previous iter
    fasta = open(out_dir(iter-1) + "output.fa", "r")
    seqs = {}
    for iter_line, line in enumerate(fasta):
        if iter_line % 2 == 0:

            if args.algorithm == 'primer':
                line_clear = re.sub(r'(_LEFT|_RIGHT)$', '', line[1:-1])

            elif args.algorithm == 'FISH':
                line_clear = line[1:-1]

            if np.isin(hash(line_clear), probe_list_hash):
                keep = True
                line_name = line[1:-1]
            else:
                keep = False
        else:
            if keep:
                seqs[line_name] = line[:-1]
    fasta.close()

    # 5*. mutations ----
    if iter != args.iterations:  # (if not final)
        if args.append:
            seqs_mutated = seqs.copy()
        else:
            seqs_mutated = dict()
        for seqs_unique in seqs.keys():
            for seqs_iter in range(args.set_size):
                init_seq = seqs[seqs_unique]
                mutated_seq = init_seq
                while init_seq == mutated_seq:
                    mutated_seq = "".join(
                        [mutate_seq(_,  mutrate=args.mutation_rate) for _ in init_seq])
                mseq = "I" + str(iter)+"N"+str(seqs_iter)+"_"
                seqs_mutated[mseq+seqs_unique] = mutated_seq

        fasta = open(out_dir(iter)+"output.fa", "w")
        for fname in seqs_mutated.keys():
            fasta.write(">" + fname + "\n" + seqs_mutated[fname]+"\n")
        fasta.close()

        merge_iter(iter)

        print("Done")

    # get merged

# </ evolutionary algorithm >
# 6. output ----

fasta = open(args.output + "/output.fa", "w")
for fname in sorted(seqs.keys()):

    if args.algorithm == 'primer':
        seqs_hits = str(probe_vals.loc[re.sub(r'(_LEFT|_RIGHT)$', '', fname)])
    elif args.algorithm == 'FISH':
        seqs_hits = str(probe_vals.loc[fname])
    fasta.write(">"+"H"+seqs_hits+"_"+fname+"\n" + seqs[fname]+"\n")

fasta.close()

write_stats(stats, args.output)
print("Done\n\nFinish")
