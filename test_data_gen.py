import argparse
import json
import sys
import random
from Bio.SeqUtils import MeltingTemp
from Bio import pairwise2
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Target and off-target sequences
    required.add_argument(
        '-nt', '--n_targets',
        type=int,
        required=True,
        help='Number of sequences in target database'
    )
    required.add_argument(
        '-no', '--n_offtargets',
        type=int,
        required=True,
        help='Number of sequences in off-target database'
    )

    optional.add_argument(
        '-minl', '--min_length',
        type=int,
        default=500,
        help='Minimal target db sequence length'
    )

    optional.add_argument(
        '-maxl', '--max_length',
        type=int,
        default=500,
        help='Maximal target db sequence length'
    )

    # Insert parameters

    optional.add_argument(
        '-il', '--insert_length',
        type=int,
        default=25,
        help='Length of an insert to generate'
    )
    optional.add_argument(
        '--insert_min_gc',
        type=float,
        default=40.0,
        help='Minimal insert GC-content'
    )
    optional.add_argument(
        '--insert_max_gc',
        type=float,
        default=60.0,
        help='Maximal insert GC-content'
    )
    optional.add_argument(
        '--insert_min_tm',
        type=float,
        default=54.0,
        help='Minimal insert Tm'
    )

    optional.add_argument(
        '--insert_max_tm',
        type=float,
        default=60.0,
        help='Maximal insert Tm'
    )

    optional.add_argument(
        '--insert_gen_attempts',
        type=int,
        default=100000,
        help='Number of attempts to generate insert with required properties.'
    )

    # Mutation rate
    optional.add_argument(
        '--insert_mutation_rate',
        type=float,
        default=0.001,
        help='Per base mutation rate of an insert when inserted into a target sequence'
    )
    optional.add_argument(
        '--offtarget_insert_mutation_rate',
        type=float,
        default=0.33,
        help='Per base mutation rate of an insert when inserted into a target sequence'
    )
    optional.add_argument(
        '--indel_prob',
        type=float,
        default=0.2,
        help='Per-base probability that a given mutation is an indel (vs SNP), range 0–1.'
    )

    # Multinsert probability
    optional.add_argument(
        '--multiinsert_prob',
        type=float,
        default=0.0,
        help='Probability of generating a second insert (0-1).'
    )

    # Technical
    optional.add_argument(
        '-t',
        '--threads',
        type=int,
        default=1,
        help='Number of threads to use.'
    )

    optional.add_argument(
        '-o', '--output_dir',
        type=Path,
        default=Path('./test_data'),
        help='Output directory for generated data.'
    )

    optional.add_argument(
        '-rs', '--random_seed',
        type=int,
        default=4,
        help='Seed for random generator.'
    )

    args = parser.parse_args()
    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)
    return args


def compute_similarity(original: str, mutated: str) -> float:
    aln = pairwise2.align.globalms(original, mutated, 1, 0, -1, -0.5, one_alignment_only=True)
    score = aln[0][2]
    return (score / len(original)) * 100


def calc_gc_content(seq: str) -> float:
    gc_count = 0
    for base in seq:
        if base == 'G' or base == 'C':
            gc_count += 1
    return (gc_count / len(seq)) * 100


def generate_insert(
        length: int,
        min_gc: float,
        max_gc: float,
        min_tm: float,
        max_tm: float,
        gen_attempts=100000) -> str:
    '''
    Generate a random DNA sequence of given length such that its GC-content (percent)
    is between min_gc and max_gc, and its melting temperature (Tm) computed via
    nearest-neighbor model is between min_tm and max_tm.
    '''

    bases = ['A', 'C', 'G', 'T']
    for _ in range(gen_attempts):  # limit attempts to avoid infinite loops
        seq = ''.join(random.choice(bases) for _ in range(length))
        gc_content = calc_gc_content(seq)
        if not (min_gc <= gc_content <= max_gc):
            continue
        # calculate Tm using nearest-neighbor model from Bio.SeqUtils
        tm = MeltingTemp.Tm_NN(seq)
        if not (min_tm <= tm <= max_tm):
            continue
        return seq
    raise ValueError(
        f"Failed to generate insert within GC[{min_gc}-{max_gc}]% and Tm[{min_tm}-{max_tm}]°C after many attempts"
    )


def mutate_seq(seq: str, mut_rate: float, indel_prob: float) -> str:
    """
    Apply mutations to a sequence: each base has `mut_rate` chance to mutate.
    With probability `indel_prob`, perform an indel (insertion or deletion), otherwise an SNP.
    """
    new_seq = []
    for base in seq:
        if random.random() < mut_rate:
            # Decide SNP vs indel
            if random.random() < indel_prob:
                # Indel: decide insertion or deletion
                if random.random() < 0.5:
                    # Insertion: add random base before original
                    new_seq.append(random.choice('ACGT'))
                    new_seq.append(base)
                else:
                    # Deletion: skip original base
                    continue
            else:
                # SNP: replace base
                choices = [b for b in 'ACGT' if b != base]
                new_seq.append(random.choice(choices))
        else:
            # No mutation: add original base
            new_seq.append(base)
    return ''.join(new_seq)


def generate_sequence(is_target: bool, insert_seq: str, args, idx: int) -> Tuple[str, str]:
    # Generate random background of variable length between min_length and max_length
    bg_len = random.randint(args.min_length, args.max_length)
    bg = ''.join(random.choice('ACGT') for _ in range(bg_len))
    seq = bg
    # Decide how many inserts (1 or possibly 2)
    n_inserts = 1 + (random.random() < args.multiinsert_prob)
    mut_rate = args.insert_mutation_rate if is_target else args.offtarget_insert_mutation_rate

    # Insert mutated copies
    print(f'Num inserts: {n_inserts}')
    for _ in range(n_inserts):
        mutated_insert = mutate_seq(insert_seq, mut_rate, args.indel_prob)
        print(mutated_insert)
        print_alignment(insert_seq, mutated_insert)
        pos = random.randint(0, len(seq))
        seq = seq[:pos] + mutated_insert + seq[pos:]
    prefix = 'target' if is_target else 'off'
    header = f">{prefix}_{idx}"
    return header, seq


def print_alignment(seq1, seq2, mode='global'):
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    alignments = aligner.align(seq1, seq2)
    print(next(alignments))

def main():
    args = parse_args()
    
if __name__ == '__main__':
    main()