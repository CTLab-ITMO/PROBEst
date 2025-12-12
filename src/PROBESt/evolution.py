import random


def mutate_position(x, mutrate, indelrate):  # issue: add indelrate to argparse
    """
    Mutates a nucleotide based on a given mutation rate and indel rate.

    This function simulates mutations in a nucleotide sequence. It can introduce
    substitutions, insertions, or deletions based on the provided mutation rate
    and indel rate.

    Args:
        x (str): The nucleotide to potentially mutate.
        indelrate (float, optional): The rate of insertions or deletions relative
                                     to the mutation rate. Defaults to 0.1.

    Returns:
        str: The mutated nucleotide or an empty string (for deletions).
    """
    # Define the set of possible nucleotides, including ambiguous codes
    nucleotide_code = "ATGCUWSMKRYBDHVN"

    # Generate a random number to determine if a mutation occurs
    rstate = random.random()

    # Check if a mutation should occur based on the mutation rate
    if rstate <= mutrate:
        # Check if the mutation should be an insertion or deletion
        if rstate <= mutrate * indelrate:
            # Insertion: Add a random nucleotide before the current one
            if rstate <= mutrate * indelrate / 2:
                return random.choice(nucleotide_code) + x
            # Deletion: Return an empty string to represent the deletion
            return ""
        # Substitution: Replace the nucleotide with a random one
        return random.choice(nucleotide_code)
    else:
        # No mutation: Return the original nucleotide
        return x


def mutate_sequence(args, out_dir, seqs):
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
                    [mutate_position(_,  mutrate=args.mutation_rate) for _ in init_seq])
            mseq = "I" + str(iter)+"N"+str(seqs_iter)+"_"
            seqs_mutated[mseq+seqs_unique] = mutated_seq

    fasta = open(out_dir+"output.fa", "w")
    for fname in seqs_mutated.keys():
        fasta.write(">" + fname + "\n" + seqs_mutated[fname]+"\n")
    fasta.close()

    print("Done")


def dedegenerate_position(x, dedegen_rate):
    """
    De-degenerates a nucleotide based on a given de-degeneration rate.

    This function attempts to replace degenerate nucleotides with specific ones.
    It's the reverse of mutation - it makes sequences less ambiguous by replacing
    degenerate codes (N, Y, W, S, M, K, R, B, D, H, V) with specific nucleotides (ATGCU).

    Args:
        x (str): The nucleotide to potentially de-degenerate.
        dedegen_rate (float): The probability of de-degenerating a degenerate nucleotide.

    Returns:
        str: The de-degenerated nucleotide or the original if not degenerate or not selected.
    """
    # Define degenerate nucleotide codes and their possible replacements
    degenerate_map = {
        'N': ['A', 'T', 'G', 'C'],
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G']
    }
    
    # If not degenerate, return as-is
    if x.upper() not in degenerate_map:
        return x
    
    # Generate a random number to determine if de-degeneration should occur
    rstate = random.random()
    
    # Check if de-degeneration should occur based on the rate
    if rstate <= dedegen_rate:
        # Replace degenerate nucleotide with a random non-degenerate option
        return random.choice(degenerate_map[x.upper()])
    
    # No de-degeneration: Return the original degenerate nucleotide
    return x
