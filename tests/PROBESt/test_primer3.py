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


import pytest
import os
import sys
from src.PROBESt.primer3 import primer_template, parse_primer3_output
from src.PROBESt.args import arguments_parse

@pytest.fixture
def mock_args(monkeypatch):
    """Mock command line arguments for testing."""
    test_args = [
        "script.py",
        "-i", "input.fasta",
        "-tb", "true_base",
        "-fb", "false_base",
        "-c", "contig_table",
        "-o", "output"
    ]
    monkeypatch.setattr(sys, 'argv', test_args)
    return arguments_parse()

@pytest.fixture
def test_dir(tmp_path):
    """Create a temporary test directory with test files."""
    # Create a sample FASTA file
    fasta_file = tmp_path / "test_sequence.exon.mod.fna"
    with open(fasta_file, "w") as f:
        f.write(">seq1\nATGCATGCATGC\n>seq2\nCGTAACGTAACG\n")

    # Create a sample Primer3 output file
    primer3_output_file = tmp_path / "primer3_output.txt"
    with open(primer3_output_file, "w") as f:
        f.write("""SEQUENCE_ID=test_sequence_seq1
PRIMER_LEFT_0_SEQUENCE=ATGC
PRIMER_RIGHT_0_SEQUENCE=CGTA
SEQUENCE_ID=test_sequence_seq2
PRIMER_LEFT_0_SEQUENCE=TTAA
PRIMER_RIGHT_0_SEQUENCE=GGCC
""")

    return {
        "fasta_file": str(fasta_file),
        "primer3_output_file": str(primer3_output_file)
    }

def test_primer_template(test_dir, mock_args):
    """Test the primer_template function."""
    # Set primer-specific arguments
    mock_args.PRIMER_PICK_PRIMER = 1
    mock_args.PRIMER_OPT_SIZE = 20
    mock_args.PRIMER_MIN_SIZE = 18
    mock_args.PRIMER_MAX_SIZE = 22
    mock_args.PRIMER_PRODUCT_SIZE_RANGE = "100-300"
    mock_args.PRIMER_NUM_RETURN = 1

    # Generate the primer template
    result = primer_template(test_dir["fasta_file"], mock_args)

    # Expected output without indentation
    expected_output = """SEQUENCE_ID=test_sequence_seq1
SEQUENCE_TEMPLATE=ATGCATGCATGC
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=100-300
PRIMER_NUM_RETURN=1
PRIMER_EXPLAIN_FLAG=1
=
SEQUENCE_ID=test_sequence_seq2
SEQUENCE_TEMPLATE=CGTAACGTAACG
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=100-300
PRIMER_NUM_RETURN=1
PRIMER_EXPLAIN_FLAG=1
="""

    # Verify the output matches the expected result
    assert result == expected_output, "Primer template output does not match expected result"

def test_parse_primer3_output(test_dir):
    """Test the parse_primer3_output function."""
    # Parse the Primer3 output file
    result = parse_primer3_output(test_dir["primer3_output_file"])

    # Expected output
    expected_output = [
        ("test_sequence_seq1", "0", "LEFT", "ATGC"),
        ("test_sequence_seq1", "0", "RIGHT", "CGTA"),
        ("test_sequence_seq2", "0", "LEFT", "TTAA"),
        ("test_sequence_seq2", "0", "RIGHT", "GGCC"),
    ]

    # Verify the parsed output matches the expected result
    assert result == expected_output, "Parsed Primer3 output does not match expected result"