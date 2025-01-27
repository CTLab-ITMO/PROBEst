import random

def mutate_seq(x, mutrate, indelrate=0.1): # issue: add indelrate to argparse
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