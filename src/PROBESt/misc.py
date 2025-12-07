import os
import re
import subprocess
import glob
from typing import List


def write_stats(stats: dict, output_dir: str):
    """
    Writes statistics about probe performance to a CSV file.

    Args:
        stats (dict): A dictionary where the keys are iteration numbers (int) and the
                      values are dictionaries containing the following keys:
                      - 'max_hits' (int): The maximum number of hits for the iteration.
                      - 'mean_hits' (float): The mean number of hits for the iteration.
        output_dir (str): The path to the output directory where the CSV file will be saved.

    Returns:
        None: The function writes the statistics to a file and does not return a value.

    Raises:
        OSError: If the output directory cannot be created or the file cannot be written.
    """
    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Construct the full path to the stats file
    stats_file = os.path.join(output_dir, 'stats.csv')

    try:
        # Write the statistics to the CSV file
        with open(stats_file, 'w') as f:
            # Write the header row
            f.write('iteration,max_hits,mean_hits\n')

            # Write each iteration's statistics
            for iter, stat in stats.items():
                f.write(f"{iter},{stat['max_hits']},{stat['mean_hits']}\n")
    except OSError as e:
        raise OSError(f"Failed to write statistics to {stats_file}: {e}")


def write_fasta(probes, output_file):
    """
    Writes primer sequences to a FASTA file.

    Args:
        probes (list): A list of tuples, where each tuple contains:
                       - sequence_id (str): The ID of the sequence.
                       - primer_num (str): The primer number.
                       - side (str): The type of primer ("LEFT" or "RIGHT").
                       - sequence (str): The primer sequence.
        output_file (str): Path to the output FASTA file.

    Returns:
        None
    """
    try:
        with open(output_file, "w") as fasta:
            for primer in probes:
                sequence_id, primer_num, side, sequence = primer
                header = f">{sequence_id}_{primer_num}_{side}"
                fasta.write(f"{header}\n{sequence}\n")
    except IOError as e:
        print(f"Error writing to file {output_file}: {e}")

def pairing(x):
    """
    Generates a pair of primer names by appending "RIGHT" and "LEFT" to a base name.

    This function takes a primer name (e.g., "primer1_LEFT" or "primer1_RIGHT") and removes
    the "LEFT" or "RIGHT" suffix to extract the base name. It then generates a pair of
    primer names by appending "RIGHT" and "LEFT" to the base name.

    Args:
        x (str): The primer name, which may contain a "LEFT" or "RIGHT" suffix.

    Returns:
        list: A list containing two strings: the base name with "RIGHT" appended and
              the base name with "LEFT" appended.
    """
    # Remove "LEFT" or "RIGHT" from the primer name
    clear_primer = re.sub(r"(RIGHT)|(LEFT)", "", string=x)

    # Return a list with the base name appended by "RIGHT" and "LEFT"
    return [clear_primer + "RIGHT", clear_primer + "LEFT"]


def collect_fasta_files(inputs: List[str]) -> List[str]:
    """
    Collect all FASTA files from input paths (files or directories).
    
    Args:
        inputs: List of input paths. Each can be:
            - A FASTA file (*.fa, *.fna, *.fasta or their .gz versions)
            - A directory containing FASTA files
    
    Returns:
        List of paths to all FASTA files found.
        
    Raises:
        ValueError: If an input path doesn't exist or no FASTA files are found.
    """
    fasta_extensions = ['.fa', '.fna', '.fasta', '.fa.gz', '.fna.gz', '.fasta.gz']
    fasta_files = []
    
    for input_path in inputs:
        if not os.path.exists(input_path):
            raise ValueError(f"Input path does not exist: {input_path}")
        
        if os.path.isfile(input_path):
            # Check if it's a FASTA file
            if any(input_path.endswith(ext) for ext in fasta_extensions):
                fasta_files.append(input_path)
            else:
                raise ValueError(f"Input file is not a recognized FASTA file: {input_path}")
        elif os.path.isdir(input_path):
            # Collect all FASTA files from directory
            dir_files = []
            for ext in fasta_extensions:
                pattern = os.path.join(input_path, f"*{ext}")
                dir_files.extend(glob.glob(pattern))
            
            if not dir_files:
                raise ValueError(f"No FASTA files found in directory: {input_path}")
            
            fasta_files.extend(sorted(dir_files))
        else:
            raise ValueError(f"Input path is neither a file nor a directory: {input_path}")
    
    if not fasta_files:
        raise ValueError("No FASTA files found in any of the provided inputs")
    
    return sorted(fasta_files)


def process_multiple_inputs(
    args: 'Namespace',
    input_fasta_files: List[str],
    output_dir: str,
    initial_generator: str = 'primer3'
) -> str:
    """
    Process multiple input FASTA files separately for initial set generation and merge their outputs.
    
    This function processes each input FASTA file independently:
    1. Creates a temporary directory for each input
    2. Runs uniline_fasta and initial_set_generation for each
    3. Merges individual outputs
    4. Combines all outputs into a single output.fa file
    
    Args:
        args: Arguments object containing all pipeline parameters
        input_fasta_files: List of paths to FASTA files to process
        output_dir: Output directory path (e.g., .tmp/0/)
        initial_generator: Initial generator to use ('primer3' or 'oligominer')
    
    Returns:
        Path to the final merged output.fa file
    
    Raises:
        ValueError: If initial_generator is unknown
        ImportError: If required modules cannot be imported
    """
    # Import here to avoid circular dependencies
    from PROBESt.bash_wrappers import uniline_fasta
    from PROBESt.merge import merge
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each input FASTA separately for initial set generation
    temp_output_dirs = []
    
    for idx, fasta_file in enumerate(input_fasta_files):
        print(f"\nProcessing input {idx + 1}/{len(input_fasta_files)}: {fasta_file}")
        
        # Create temporary directory for this input
        temp_dir = os.path.join(output_dir, f"input_{idx}/")
        os.makedirs(temp_dir, exist_ok=True)
        temp_output_dirs.append(temp_dir)
        
        # Create a temporary args object with this specific input file
        temp_args = type('Args', (), {})()
        for attr in dir(args):
            if not attr.startswith('_'):
                setattr(temp_args, attr, getattr(args, attr))
        temp_args.input = fasta_file
        
        # Make uniline fasta for this input
        uniline_fasta(temp_args, temp_dir)
        print(f"Input fasta parsed: {fasta_file}")
        
        # Template generation for this input
        if initial_generator == "primer3":
            from PROBESt.primer3 import initial_set_generation
            initial_set_generation(temp_args, temp_dir)
        elif initial_generator == "oligominer":
            from PROBESt.oligominer import initial_set_generation
            initial_set_generation(temp_args, temp_dir)
        else:
            raise ValueError(f"Unknown initial generator: {initial_generator}")
        
        # Merge this input's output
        merge(algo=temp_args.algorithm,
              input=os.path.join(temp_dir, "output.fa"),
              output=os.path.join(temp_dir, "merged.fa"),
              tmp=os.path.join(temp_dir, "fasta_table.tsv"),
              NNN=10,
              script_path=temp_args.script_path)
    
    # Merge all outputs to the first directory's output.fa
    print("\nMerging all initial set outputs...")
    final_output_fa = os.path.join(output_dir, "output.fa")
    with open(final_output_fa, 'w') as outfile:
        for temp_dir in temp_output_dirs:
            merged_fa = os.path.join(temp_dir, "merged.fa")
            if os.path.exists(merged_fa) and os.path.getsize(merged_fa) > 0:
                with open(merged_fa, 'r') as infile:
                    outfile.write(infile.read())
            else:
                # If merged.fa doesn't exist, try output.fa
                output_fa = os.path.join(temp_dir, "output.fa")
                if os.path.exists(output_fa) and os.path.getsize(output_fa) > 0:
                    with open(output_fa, 'r') as infile:
                        outfile.write(infile.read())
    
    return final_output_fa
