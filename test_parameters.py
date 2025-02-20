'''
Script for running the main pipeline with a given parameters grid to collect hit statistics.
'''

import argparse
import itertools
import json
import os
import pandas as pd
import subprocess
import tqdm

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from math import prod
from sys import argv, exit


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        '-p',
        '--params_grid',
        type=Path,
        required=True,
        help='JSON file with a parameters grid for benchmark.'
    )

    optional.add_argument(
        '-o',
        '--output',
        type=Path,
        default='./data/param_search',
        help='Output files directory.'
    )

    optional.add_argument(
        '-t',
        '--threads',
        type=int,
        default=1,
        help='Number of threads'
    )

    optional.add_argument(
        '-tpw',
        '--threads_per_worker',
        type=int,
        default=1,
        help='Number of threads per pipeline.'
    )

    optional.add_argument(
        '-kf',
        '--keep_failed',
        action='store_true',
        help='Include failed runs in the final table'
    )

    args = parser.parse_args()
    if len(argv) < 2:
        parser.print_usage()
        exit(1)
    return args


def read_param_grid(param_file) -> dict:
    with open(param_file) as json_grid:
        return json.load(json_grid)


def create_if_not_exist(dir_path: Path) -> None:
    '''Creates provided directory if it doesn't exist'''
    if not dir_path.exists():
        os.mkdir(dir_path)


def read_stats(output_dir: str) -> dict:
    '''Reads the stats and returns dict with max and mean hits for the first and last iteration'''

    stats_file = Path(output_dir) / 'stats.csv'
    stats = pd.read_csv(stats_file)

    return stats


def run_pipeline(params: dict, run_number: int) -> list:

    # output - the only obligatory key in the params dict
    pipeline_file = Path(__file__).parent / 'pipeline.py'
    output_dir = params['output']
    create_if_not_exist(output_dir)

    # Construct run command
    command = f'python {pipeline_file} '
    for param, value in params.items():
        if isinstance(value, list):
            value = ' '.join(map(str, value))
        command += f'--{param} {value} '

    # Run pipeline
    log_file = Path(output_dir) / 'pipeline.log'
    with open(log_file, 'w') as log:
        subprocess.run(command, stdout=log, shell=True, executable="/bin/bash")

    stats = read_stats(output_dir)
    output = []
    for _, row in stats.iterrows():
        iteration_data = {
            'run': run_number,
            'iteration': row['iteration'],
            'max_hits': row['max_hits'],
            'mean_hits': row['mean_hits'],
        }
        iteration_data.update(params)
        output.append(iteration_data)
    return output


def write_output(param_stats: list, output_dir: str) -> None:
    output_file = Path(output_dir) / 'param_stats.csv'

    df = pd.DataFrame(param_stats).sort_values(by=['run', 'iteration'])

    # Convert lists to strings
    for column in df.columns:
        df[column] = df[column].apply(lambda x: ' '.join(map(str, x)) if isinstance(x, list) else x)

    df.to_csv(output_file, index=False)


def main() -> None:

    # Processing input
    args = parse_args()
    params_grid = read_param_grid(args.params_grid)
    param_combinations = itertools.product(*params_grid.values())  # Get all possible combinations of parameters
    comb_num = prod(len(p) for p in params_grid.values())  # Calculate number of combinations

    # Format parameters for parallel execution
    params_for_executor = []
    for c, combination in enumerate(param_combinations):
        output_directory = args.output / f'{c}/'

        params_values = list(combination)
        params_values.append(output_directory)
        param_names = list(params_grid.keys())
        param_names.append('output')
        params_final = {
            param: value
            for param, value in zip(param_names, params_values)
        }
        params_for_executor.append(params_final)

    # Test params in parallel
    create_if_not_exist(args.output)
    results = []
    num_workers = args.threads // args.threads_per_worker

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        run_counter = itertools.count(start=1)
        futures = {}

        for params in params_for_executor:
            run_number = next(run_counter)
            future = executor.submit(run_pipeline, params, run_number)
            futures[future] = params

        for future in tqdm.tqdm(as_completed(futures), total=comb_num, desc='Progress'):
            params = futures[future]
            try:
                result = future.result()
                results.extend(result)
            except Exception as e:
                na_stats = {
                    "iteration": 1,
                    "max_hits": 0,
                    "mean_hits": 0
                }
                if args.keep_failed:
                    results.append({**params, **na_stats})
                print(f'Error with params: {params}: {e}')

    write_output(results, args.output)


if __name__ == '__main__':
    main()
