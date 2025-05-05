import argparse
import sys
import random
from pathlib import Path
from Bio.SeqUtils import MeltingTemp


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
        '-ml', '--min_length',
        type=int,
        default=500,
        help='Minimal sequence length'
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

if __name__ == '__main__':
    args = parse_args()
    print(generate_insert(args.insert_length, args.insert_min_gc, args.insert_max_gc, args.insert_min_tm, args.insert_max_tm, gen_attempts=args.insert_gen_attempts))