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


from Bio import SeqIO
import os
import sys


def primer_template(fasta_file: str) -> str:
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
        seq_id = os.path.basename(fasta_file).replace(".fasta", "")

        # Generate the primer design template for the current sequence
        template = f"""SEQUENCE_ID={seq_id}
SEQUENCE_TEMPLATE={record.seq}
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=10
PRIMER_PICK_RIGHT_PRIMER=10
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE=25
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=30
PRIMER_PRODUCT_SIZE_RANGE=100-1000
PRIMER_NUM_RETURN=10
PRIMER_EXPLAIN_FLAG=1
="""
        output.append(template)

    # Join all templates into a single string and return
    return "\n".join(output)


if __name__ == '__main__':
    fasta = sys.argv[1]
    primer_temp = primer_template(fasta)
    template_path = './primer3_template'
    with open(template_path, "w") as template:
        template.writelines(primer_temp)