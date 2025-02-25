from Bio import SeqIO
import os
import subprocess


def primer_template(fasta_file: str, args) -> str:
    """
    Generates a primer design template for each sequence in a given FASTA file.

    This function reads sequences from a FASTA file and generates a primer design template
    for each sequence. The template is formatted for use with primer design tools like Primer3.

    Args:
        fasta_file (str): Path to the input FASTA file containing sequences.
        args: Parsed command-line arguments containing Primer3 template options.

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
PRIMER_PICK_LEFT_PRIMER={args.PRIMER_PICK_PRIMER}
PRIMER_PICK_RIGHT_PRIMER={args.PRIMER_PICK_PRIMER}
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE={args.PRIMER_OPT_SIZE}
PRIMER_MIN_SIZE={args.PRIMER_MIN_SIZE}
PRIMER_MAX_SIZE={args.PRIMER_MAX_SIZE}
PRIMER_PRODUCT_SIZE_RANGE={args.PRIMER_PRODUCT_SIZE_RANGE}
PRIMER_NUM_RETURN={args.PRIMER_NUM_RETURN}
PRIMER_EXPLAIN_FLAG=1
="""
        output.append(template)

    # Join all templates into a single string and return
    return "\n".join(output)


def parse_primer3_output(primer3_file: str) -> list:
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


def primer2fasta(args, out_dir: str, probes: list) -> None:
    """
    Writes primer sequences to a FASTA file and optionally appends additional sequences.

    This function writes the extracted primer sequences to a FASTA file. If additional
    sequences are provided via `args.add_set`, they are appended to the FASTA file.

    Args:
        args: Parsed command-line arguments.
        out_dir (str): Path to the output directory.
        probes (list): A list of tuples containing primer information (sequence_id, primer_num, side, sequence).
    """
    # Write primers to the output FASTA file
    output_fasta_path = os.path.join(out_dir, "output.fa")
    with open(output_fasta_path, "w") as fasta:
        for probe in probes:
            sequence_id, primer_num, side, sequence = probe
            header = f">{sequence_id}_{primer_num}_{side}"
            fasta.write(f"{header}\n{sequence}\n")

    # Append additional sequences if provided
    if args.add_set:
        add_fasta_cmd = f"cat {args.add_set} >> {output_fasta_path}"
        subprocess.run(add_fasta_cmd, shell=True, executable="/bin/bash")


def initial_set_generation(args, out_dir: str) -> None:
    """
    Generates an initial set of primers using Primer3 and writes them to a FASTA file.

    This function performs the following steps:
    1. Generates a Primer3 template from the input FASTA file.
    2. Runs Primer3 to generate primer sequences.
    3. Parses the Primer3 output and writes the primers to a FASTA file.
    4. Optionally appends additional sequences to the FASTA file.

    Args:
        args: Parsed command-line arguments.
        out_dir (str): Path to the output directory.
    """
    # Generate the Primer3 template
    primer_temp = primer_template(os.path.join(out_dir, "input.fa"), args)

    # Write the template to a file
    template_path = os.path.join(out_dir, "template")
    with open(template_path, "w") as template:
        template.writelines(primer_temp)

    # Run Primer3
    primer3_output_path = os.path.join(out_dir, "output.p3")
    primer3_cmd = f"{args.primer3} {template_path} --output {primer3_output_path}"
    subprocess.run(primer3_cmd, shell=True, executable="/bin/bash")
    print("Primer3 done")

    # Parse the Primer3 output and write to FASTA
    probes = parse_primer3_output(primer3_output_path)
    primer2fasta(args, out_dir, probes)
