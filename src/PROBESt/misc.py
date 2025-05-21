import os
import re
import subprocess


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
