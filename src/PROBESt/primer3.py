from Bio import SeqIO
import os

def primer_template(fasta_file,
                    PRIMER_PICK_PRIMER,
                    PRIMER_OPT_SIZE,
                    PRIMER_MIN_SIZE,
                    PRIMER_MAX_SIZE,
                    PRIMER_PRODUCT_SIZE_RANGE,
                    PRIMER_NUM_RETURN):
    """
    Generates a primer design template for each sequence in a given FASTA file.

    This function reads sequences from a FASTA file and generates a primer design template
    for each sequence. The template is formatted for use with primer design tools like Primer3.

    Args:
        fasta_file (str): Path to the input FASTA file containing sequences.
        PRIMER_PICK_PRIMER (int): Flag to indicate whether to pick primers (1 for yes, 0 for no).
        PRIMER_OPT_SIZE (int): Optimal primer length.
        PRIMER_MIN_SIZE (int): Minimum primer length.
        PRIMER_MAX_SIZE (int): Maximum primer length.
        PRIMER_PRODUCT_SIZE_RANGE (str): Desired size range for the PCR product (e.g., "100-300").
        PRIMER_NUM_RETURN (int): Number of primer pairs to return.

    Returns:
        str: A concatenated string of primer design templates for all sequences in the FASTA file.
    """
    
    # Read sequences from the FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    output = []

    # Iterate over each sequence in the FASTA file
    for record in sequences:
        # Extract the base filename without the extension and modifier
        seq_id = os.path.basename(fasta_file).replace(".exon.mod.fna", "")

        # Generate the primer design template for the current sequence
        template = f"""SEQUENCE_ID={seq_id}_{record.id}
SEQUENCE_TEMPLATE={record.seq}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER={PRIMER_PICK_PRIMER}
PRIMER_PICK_RIGHT_PRIMER={PRIMER_PICK_PRIMER}
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE={PRIMER_OPT_SIZE}
PRIMER_MIN_SIZE={PRIMER_MIN_SIZE}
PRIMER_MAX_SIZE={PRIMER_MAX_SIZE}
PRIMER_PRODUCT_SIZE_RANGE={PRIMER_PRODUCT_SIZE_RANGE}
PRIMER_NUM_RETURN={PRIMER_NUM_RETURN}
PRIMER_EXPLAIN_FLAG=1
="""
        output.append(template)
    
    # Join all templates into a single string and return
    return "\n".join(output)

def parse_primer3_output(primer3_file):
    """
    Parses a Primer3 output file and extracts primer sequences along with their metadata.

    This function reads a Primer3 output file and extracts information about the primers,
    including their sequence ID, primer number, type (LEFT or RIGHT), and sequence.

    Args:
        primer3_file (str): Path to the Primer3 output file.

    Returns:
        list: A list of tuples, where each tuple contains:
              - sequence_id (str): The ID of the sequence.
              - primer_num (str): The primer number (e.g., "0" for the first primer pair).
              - primer_type (str): The type of primer ("LEFT" or "RIGHT").
              - primer_seq (str): The primer sequence.
    """
    probes = []
    sequence_id = None

    with open(primer3_file, "r") as file:
        for line in file:
            line = line.strip()

            # Extract SEQUENCE_ID
            if line.startswith("SEQUENCE_ID="):
                sequence_id = line.split("=")[1]

            # Extract LEFT primer sequence
            if line.startswith("PRIMER_LEFT_") and "SEQUENCE" in line:
                primer_num = line.split("_")[2]
                primer_seq = line.split("=")[1]
                probes.append((sequence_id, primer_num, "LEFT", primer_seq))

            # Extract RIGHT primer sequence
            if line.startswith("PRIMER_RIGHT_") and "SEQUENCE" in line:
                primer_num = line.split("_")[2]
                primer_seq = line.split("=")[1]
                probes.append((sequence_id, primer_num, "RIGHT", primer_seq))

    return probes
