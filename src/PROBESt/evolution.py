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
