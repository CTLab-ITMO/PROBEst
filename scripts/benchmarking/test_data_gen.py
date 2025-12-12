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
import json
import pandas as pd
import random
import sys
from Bio.SeqUtils import MeltingTemp
from Bio.Align import PairwiseAligner
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple, List


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
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(original, mutated)
    score = alignments[0].score
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
        if not (min_gc <= calc_gc_content(seq) <= max_gc):
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


def generate_sequence(
        is_target: bool, insert_seq: str, args, idx: int, seed: int
        ) -> Tuple[Tuple[str, str], List[Tuple]]:

    random.seed(seed)
    # Generate random background of variable length between min_length and max_length
    bg_len = random.randint(args.min_length, args.max_length)
    bg = ''.join(random.choice('ACGT') for _ in range(bg_len))
    seq = bg
    # Decide how many inserts (1 or possibly 2)
    n_inserts = 1 + (random.random() < args.multiinsert_prob)
    mut_rate = args.insert_mutation_rate if is_target else args.offtarget_insert_mutation_rate
    metadata = []

    # Insert mutated copies
    for ins_idx in range(n_inserts):
        mutated = mutate_seq(insert_seq, mut_rate, args.indel_prob)
        pos = random.randint(0, len(seq))
        seq = seq[:pos] + mutated + seq[pos:]
        rec_id = f"{'target' if is_target else 'off'}_{idx}"
        sim = compute_similarity(insert_seq, mutated)
        metadata.append({
            'record_id': rec_id,
            'insert_idx': ins_idx + 1,
            'inserted_seq': mutated,
            'position': pos,
            'similarity': sim
        })
    header = f">{'target' if is_target else 'off'}_{idx}"
    return (header, seq), metadata


def main():
    args = parse_args()
    random.seed(args.random_seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # write all params
    with open(args.output_dir / 'params.json', 'w') as pf:
        json.dump(
            {
                k: v for k, v in vars(args).items()
                if (k != 'metadata' and not isinstance(v, Path))
            }, pf, indent=2)


    print('Generating insert')
    # generate initial insert and save it's parameters
    insert_seq = generate_insert(
        args.insert_length,
        args.insert_min_gc,
        args.insert_max_gc,
        args.insert_min_tm,
        args.insert_max_tm)

    insert_info = {
        'sequence': insert_seq,
        'length': len(insert_seq),
        'gc_content': calc_gc_content(insert_seq),
        'tm': MeltingTemp.Tm_NN(insert_seq)
    }

    with open(args.output_dir / 'insert.json', 'w') as inf:
        json.dump(insert_info, inf, indent=2)

    # parallel generation
    print('Generating databases')
    records: List[Tuple[str, str]] = []
    metadata: List[dict] = []
    base_seed = args.random_seed
    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        futures = []
        # target sequences
        for i in range(args.n_targets):
            futures.append(pool.submit(generate_sequence, True, insert_seq, args, i, base_seed + i))
        # off-target sequences
        for j in range(args.n_offtargets):
            futures.append(pool.submit(generate_sequence, False, insert_seq, args, j, base_seed + 10000 + j))

    for fut in as_completed(futures):
        (header, seq), meta = fut.result()
        records.append((header, seq))
        metadata.extend(meta)

    print('Writing output')
    # Write FASTA
    tgt_dir = args.output_dir / 'target'
    off_dir = args.output_dir / 'offtarget'
    tgt_dir.mkdir(exist_ok=True)
    off_dir.mkdir(exist_ok=True)
    for header, seq in records:
        rid = header.lstrip('>')
        subdir = tgt_dir if rid.startswith('target_') else off_dir
        with open(subdir/f"{rid}.fasta", 'w') as f:
            f.write(f"{header}\n{seq}\n")

    # Create DataFrame from metadata and write TSV
    df = pd.DataFrame(metadata)
    df.to_csv(args.output_dir/'metadata.tsv', sep='\t', index=False)

    # Compute metrics using pandas
    # Count sequences with multi-inserts
    target_insert_counts = df[df['record_id'].str.startswith('target_')].groupby('record_id').size()
    off_insert_counts = df[df['record_id'].str.startswith('off_')].groupby('record_id').size()
    multi_target = (target_insert_counts > 1).sum()
    multi_off = (off_insert_counts > 1).sum()
    # Average similarity
    avg_target = df[df['record_id'].str.startswith('target_')]['similarity'].mean()
    avg_off = df[df['record_id'].str.startswith('off_')]['similarity'].mean()

    print('Finished!')
    # Print summary
    print('Summary:')
    print(f"Number of target sequences with multi-inserts: {multi_target}")
    print(f"Number of off-target sequences with multi-inserts: {multi_off}")
    print(f"Average similarity in target_db: {avg_target:.2f}%")
    print(f"Average similarity in off_target_db: {avg_off:.2f}%")


if __name__ == '__main__':
    main()
