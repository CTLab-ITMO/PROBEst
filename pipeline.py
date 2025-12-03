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
#
# < de-degeneration algorithm >
# 7. de-degeneration iterations
# </ de-degeneration algorithm >
#
# 8. modeling and visualization

# 0. Imports: python packages ----

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import re

# Add src directory to path to ensure local modules can be imported
script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(script_dir, 'src')
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

# 0 Imports: absolute import ----
from PROBESt.primer3 import initial_set_generation
from PROBESt.evolution import mutate_position
from PROBESt.bash_wrappers import uniline_fasta, blastn_function, probe_check_function
from PROBESt.misc import write_stats
from PROBESt.merge import merge
from PROBESt.args import arguments_parse
from PROBESt.modeling import run_modeling
from PROBESt.dedegeneration import run_dedegeneration

# Functions


def out_dir(iter: int):
    return args.output_tmp + str(iter) + "/"


def merge_iter(iter: int):
    out = out_dir(iter)
    merge(algo=args.algorithm,
          input=out + "output.fa",
          output=out + "merged.fa",
          tmp=out + "fasta_table.tsv",
          NNN=10,
          script_path=args.script_path)


# 0. Argparsing ----
args = arguments_parse()
args.script_path = os.path.dirname(
    os.path.realpath(__file__)) + \
    "/scripts/generator/"

# 1. Initial set generation ----
print("\n---- PROBESt v.0.2.0 ----\n")
print("Arguments passed")

# Create TMP
os.makedirs(out_dir(0), exist_ok=True)

# Make uniline fasta
uniline_fasta(args, out_dir(0))
print("Input fasta parsed")

# Template generation
initial_set_generation(args, out_dir(0))
merge_iter(0)

# < evolutionary algorithm >

# commands
blastn = blastn_function(args)
probe_check = probe_check_function(args)

# arrays
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
                        [mutate_position(_,  mutrate=args.mutation_rate, indelrate=args.indel_rate) for _ in init_seq])
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

# 7. De-degeneration ----
if args.dedegeneration_iterations > 0:
    dedegeneration_input = args.output + "/output.fa"
    dedegeneration_output = args.output + "/output_dedegenerated.fa"
    run_dedegeneration(args, dedegeneration_input, dedegeneration_output)
    # Use de-degenerated output for modeling
    final_output_fa = dedegeneration_output
else:
    final_output_fa = args.output + "/output.fa"

# 8. Modeling and visualization ----
modeling_output = os.path.join(args.output, "modeling_results.tsv")
run_modeling(args, args.input, final_output_fa, modeling_output)
