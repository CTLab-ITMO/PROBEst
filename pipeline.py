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
from PROBESt.evolution import mutate_position
from PROBESt.bash_wrappers import uniline_fasta, blastn_function, probe_check_function
from PROBESt.misc import write_stats
from PROBESt.merge import merge
from PROBESt.args import arguments_parse
from PROBESt.modeling import run_modeling
from PROBESt.dedegeneration import run_dedegeneration
from PROBESt.prepare_blast import prepare_bases_if_needed

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

# Create output directory early (needed for BLAST database preparation)
os.makedirs(args.output, exist_ok=True)

# 0.1. Prepare BLAST databases from FASTA directories if needed ----
args.true_base, args.false_base = prepare_bases_if_needed(
    true_base=args.true_base,
    false_bases=args.false_base,
    output_dir=args.output,
    contig_table_path=args.contig_table,
    tmp_dir=args.prep_db_tmp,
    script_path=args.script_path
)

# Create TMP
os.makedirs(out_dir(0), exist_ok=True)

# Make uniline fasta
uniline_fasta(args, out_dir(0))
print("Input fasta parsed")

# Template generation
# Check if initial_generator argument exists (for backward compatibility)
initial_generator = getattr(args, 'initial_generator', 'primer3')

if initial_generator == "primer3":
    from PROBESt.primer3 import initial_set_generation
    initial_set_generation(args, out_dir(0))
elif initial_generator == "oligominer":
    from PROBESt.oligominer import initial_set_generation
    initial_set_generation(args, out_dir(0))
else:
    raise ValueError(f"Unknown initial generator: {initial_generator}")
merge_iter(0)

# Check if initial set generation produced any probes
output_fa_path = out_dir(0) + "output.fa"
if not os.path.exists(output_fa_path) or os.path.getsize(output_fa_path) == 0:
    print("\nERROR: Initial set generation produced no probes (empty output.fa).")
    print("This may happen if:")
    print("  - The input sequences are too short or have low complexity")
    print("  - The probe generation parameters are too restrictive")
    print("  - The selected initial generator is not suitable for your input")
    print(f"\nSuggestion: Try using a different initial generator (current: {initial_generator})")
    if initial_generator == "primer3":
        print("  Consider using: --initial_generator oligominer")
    elif initial_generator == "oligominer":
        print("  Consider using: --initial_generator primer3")
    print("\nAlternatively, adjust the probe generation parameters:")
    if initial_generator == "oligominer":
        print("  - Try adjusting --oligominer_probe_length")
        print("  - Try adjusting --oligominer_temperature")
    sys.exit(1)

# < evolutionary algorithm >

# commands
blastn = blastn_function(args)
probe_check = probe_check_function(args)

# arrays
stats = {}  # Hit stats for each iteration
last_valid_iter = 0  # Track the last iteration with valid probes

for iter in range(1, args.iterations+1):
    print("\nIteration", iter, "----")
    os.makedirs(out_dir(iter), exist_ok=True)

    # Check if previous iteration's merged.fa exists and is not empty
    prev_merged_fa = out_dir(iter-1) + "merged.fa"
    if not os.path.exists(prev_merged_fa) or os.path.getsize(prev_merged_fa) == 0:
        if iter == 1:
            # This means iteration 0 had empty output, which should have been caught earlier
            # But handle it here as a safety check
            print("\nERROR: No probes available from previous iteration.")
            print("Initial set generation may have failed. Try using a different initial generator.")
            sys.exit(1)
        else:
            # Later iteration had empty output - use last valid iteration
            print(f"\nWARNING: Iteration {iter-1} produced no probes (empty merged.fa).")
            print(f"Stopping evolutionary algorithm. Using results from iteration {last_valid_iter}.")
            break

    # 2. blastn ----
    blastn_iter = blastn + " -query " + prev_merged_fa

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
    clear_hits_path = out_dir(iter) + "clear_hits.tsv"
    try:
        # Check if file exists and is not empty
        if not os.path.exists(clear_hits_path) or os.path.getsize(clear_hits_path) == 0:
            if iter == 1:
                # First iteration failed - this is a problem with initial set or probe_check
                raise InterruptedError(
                    "Empty file after filtration at iteration 1. "
                    "This suggests the initial probe set may not be suitable. "
                    "Try using a different initial generator (--initial_generator) or "
                    "adjusting probe_check properties and reviewing false databases.")
            else:
                # Later iteration failed - use last valid iteration
                print(f"\nWARNING: Iteration {iter} produced no valid probes after filtration.")
                print(f"Stopping evolutionary algorithm. Using results from iteration {last_valid_iter}.")
                break
        
        probe_out = pd.read_table(clear_hits_path, sep=' ', header=None)
    except (pd.errors.EmptyDataError, FileNotFoundError) as e:
        if iter == 1:
            raise InterruptedError(
                "Empty file after filtration at iteration 1. "
                "This suggests the initial probe set may not be suitable. "
                "Try using a different initial generator (--initial_generator) or "
                "adjusting probe_check properties and reviewing false databases.")
        else:
            print(f"\nWARNING: Iteration {iter} produced no valid probes after filtration.")
            print(f"Stopping evolutionary algorithm. Using results from iteration {last_valid_iter}.")
            break

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
    prev_output_fa = out_dir(iter-1) + "output.fa"
    if not os.path.exists(prev_output_fa) or os.path.getsize(prev_output_fa) == 0:
        if iter == 1:
            print("\nERROR: No probes available from previous iteration.")
            print("Initial set generation may have failed. Try using a different initial generator.")
            sys.exit(1)
        else:
            print(f"\nWARNING: Iteration {iter-1} produced no probes (empty output.fa).")
            print(f"Stopping evolutionary algorithm. Using results from iteration {last_valid_iter}.")
            break
    
    fasta = open(prev_output_fa, "r")
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
    
    # Check if we have any sequences after filtering
    if len(seqs) == 0:
        if iter == 1:
            print("\nERROR: No probes matched after filtering at iteration 1.")
            print("This suggests the initial probe set may not be suitable.")
            print("Try using a different initial generator (--initial_generator).")
            sys.exit(1)
        else:
            print(f"\nWARNING: Iteration {iter} produced no matching probes after filtering.")
            print(f"Stopping evolutionary algorithm. Using results from iteration {last_valid_iter}.")
            break
    
    # Mark this iteration as valid
    last_valid_iter = iter

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
        
        # Check if merged.fa was created and is not empty
        merged_fa_path = out_dir(iter) + "merged.fa"
        if not os.path.exists(merged_fa_path) or os.path.getsize(merged_fa_path) == 0:
            print(f"\nWARNING: Iteration {iter} produced empty merged.fa after mutations.")
            print(f"Stopping evolutionary algorithm. Using results from iteration {last_valid_iter}.")
            break

    print("Done")

    # get merged

# </ evolutionary algorithm >
# 6. output ----
# Use the last valid iteration's data for final output
if last_valid_iter == 0:
    # No valid iterations - this should have been caught earlier, but handle it
    print("\nERROR: No valid iterations completed.")
    sys.exit(1)

# If we broke out early, we need to reload the data from the last valid iteration
# Note: last_valid_iter is the iteration that had valid probes after filtering
# We need to use the probe_vals and seqs from that iteration
final_iter = last_valid_iter
if final_iter < args.iterations:
    # We broke out early, reload data from the last valid iteration
    print(f"\nReloading data from iteration {final_iter} (last valid iteration).")
    clear_hits_path = out_dir(final_iter) + "clear_hits.tsv"
    probe_out = pd.read_table(clear_hits_path, sep=' ', header=None)
    probe_vals = probe_out.iloc[:, 0].value_counts()
    probe_list = list(set(probe_out.iloc[:, 0]))
    probe_list_hash = [hash(_) for _ in probe_list]
    
    # Reload sequences from the iteration before the last valid one
    # (the output.fa that was used to generate the last valid iteration's results)
    prev_output_fa = out_dir(final_iter-1) + "output.fa"
    fasta = open(prev_output_fa, "r")
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
    
    print(f"Using {len(seqs)} probes from iteration {final_iter-1} filtered by iteration {final_iter}.")
else:
    # We completed all iterations - seqs and probe_vals are already set from the last iteration
    print(f"\nCompleted all {args.iterations} iterations. Using final results.")

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
