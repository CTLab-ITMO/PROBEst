import os
import re

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
        
def out_dir(iter, output, output_tmp):
    """
    Warning: deprecated
        
    Generates a directory path for temporary or output files based on the iteration number.

    Args:
        iter (int): The iteration number, used to create a unique subdirectory.

    Returns:
        str: The generated directory path.
    """
    try:
        if args.output_tmp == "":
            # Use the default output directory with a `.tmp` subdirectory
            path = os.path.join(args.output, ".tmp", str(iter))
        else:
            # Use the custom temporary output directory
            path = os.path.join(args.output_tmp, str(iter))

        # Create the directory if it does not exist
        os.makedirs(path, exist_ok=True)
        return path + "/"
    except AttributeError:
        raise AttributeError("The `args` object is missing required attributes (`output` or `output_tmp`).")
    
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