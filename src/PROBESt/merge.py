import subprocess
from shutil import copyfile

def merge(algo, input, output, tmp, NNN, script_path):
    '''
    Merges primer sequences from a FASTA file based on the specified algorithm.

    Parameters:
    algo (str): The algorithm to use for merging; parse form args.algorithm.
    input (str): Path to the input FASTA file containing primer pairs. 
                 Pairs should be specified as 'LEFT' and 'RIGHT'.
    output (str): Path to the output file where merged sequences will be written.
    tmp (str): Path to a temporary file for storing the converted FASTA table.
    NNN (int): The number of 'N' characters to use as a separator between sequences.
    script_path (str): The directory path where bash scripts are located.

    Raises:
    LookupError: If there is an inconsistency in sequence names in the FASTA file.
    KeyError: If an unsupported algorithm is specified.
    '''

    # Construct the command to convert FASTA to a table format using a bash script
    fasta2table_command = f"bash {script_path}fasta2table.sh -i {input} -o {tmp}"

    # Execute the bash command
    subprocess.run(fasta2table_command, shell=True)

    if algo == "primer":
        # Open the temporary FASTA table for reading
        with open(tmp, "r") as fasta_table, open(output, "w") as output_fasta:
            sepN = "N" * NNN  # Create a separator string of 'N's

            # Iterate through lines in the FASTA table
            for i, line in enumerate(fasta_table):
                dev = i % 4  # Determine the line type based on its position
                if dev == 0:
                    inp_name = line[:-6]  # Extract input name (remove trailing characters)
                elif dev == 1:
                    left = line.strip()  # Get the left primer sequence (remove newline)
                elif dev == 2:
                    # Verify that sequence names match; raise an error if not
                    if hash(inp_name) != hash(line[0:(len(line)-7)]):
                        raise LookupError(f"Error on line {i}: check names in FASTA file")
                        break
                else:
                    right = line.strip()  # Get the right primer sequence
                    # Write the merged sequence to the output file
                    output_fasta.write(f"{inp_name}\n{left}{sepN}{right}\n")

        print("Primers merged successfully.")

    elif algo == "FISH":
        # For FISH algorithm, simply copy the temporary file to output
        copyfile(tmp, output)

    else:
        # Raise an error if an unsupported algorithm is provided
        raise KeyError("Algorithm not implemented yet; see docs or help.")

