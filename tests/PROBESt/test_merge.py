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
import subprocess
from src.PROBESt.merge import merge

@pytest.fixture
def mock_script(tmp_path):
    """Create a mock fasta2table.sh script for testing."""
    script_dir = tmp_path / "scripts"
    os.makedirs(script_dir)
    script_path = script_dir / "fasta2table.sh"
    
    # Create a simple mock script that copies input to output
    with open(script_path, "w") as f:
        f.write("""#!/bin/bash
while getopts "i:o:" opt; do
  case $opt in
    i) input="$OPTARG";;
    o) output="$OPTARG";;
  esac
done

# Simple mock implementation: just copy input to output
cp "$input" "$output"
""")
    
    # Make the script executable
    os.chmod(script_path, 0o755)
    return str(script_dir) + "/"

@pytest.fixture
def test_dir(tmp_path):
    """Create a temporary test directory with test files."""
    # Create a sample input FASTA file for primer merging
    input_fasta_primer = tmp_path / "input_primer.fasta"
    with open(input_fasta_primer, "w") as f:
        f.write(">seq1_LEFT\nATGC\n>seq1_RIGHT\nCGTA\n>seq2_LEFT\nTTAA\n>seq2_RIGHT\nGGCC\n")
    input_fasta_primer3 = tmp_path / "input_primer3_style.fasta"
    with open(input_fasta_primer3, "w") as f:
        # Same pattern as primer3.primer2fasta: {id}_{pair}_{LEFT|RIGHT}
        f.write(
            ">contigA_0_LEFT\nATGC\n>contigA_0_RIGHT\nCGTA\n"
            ">contigA_1_LEFT\nTTAA\n>contigA_1_RIGHT\nGGCC\n"
        )

    # Create a sample input FASTA file for FISH merging
    input_fasta_fish = tmp_path / "input_fish.fasta"
    with open(input_fasta_fish, "w") as f:
        f.write(">seq1\nATGC\n>seq2\nCGTA\n")

    # Define output and temporary files
    output_primer = tmp_path / "output_primer.fasta"
    output_fish = tmp_path / "output_fish.fasta"
    tmp_file = tmp_path / "tmp_table.txt"

    return {
        "input_fasta_primer": str(input_fasta_primer),
        "input_fasta_primer3": str(input_fasta_primer3),
        "input_fasta_fish": str(input_fasta_fish),
        "output_primer": str(output_primer),
        "output_fish": str(output_fish),
        "tmp_file": str(tmp_file),
    }

def test_merge_primer_primer3_indexed_headers(test_dir, mock_script):
    """Primer3 FASTA uses {id}_{pair}_LEFT / _RIGHT (primer3.primer2fasta)."""
    out = str(test_dir["output_primer"]) + "_p3"
    tmp = str(test_dir["tmp_file"]) + "_p3"
    merge("primer", test_dir["input_fasta_primer3"], out, tmp, 5, mock_script)
    with open(out, "r") as f:
        body = f.read()
    assert body == ">contigA_0\nATGCNNNNNCGTA\n>contigA_1\nTTAANNNNNGGCC\n"


def test_merge_primer_already_merged_is_passthrough(test_dir, mock_script, tmp_path):
    """Evolution / dedeg output: one record per probe — copy without pair merge."""
    inp = tmp_path / "single.fa"
    inp.write_text(">probe1\nAAAAANNNNNBBBB\n", encoding="utf-8")
    out = tmp_path / "out.fa"
    merge("primer", str(inp), str(out), str(tmp_path / "t.tsv"), 5, mock_script)
    assert out.read_text(encoding="utf-8") == ">probe1\nAAAAANNNNNBBBB\n"


def test_merge_primer_algorithm(test_dir, mock_script):
    """Test the merge function with the 'primer' algorithm."""
    # Define input parameters
    algo = "primer"
    input_file = test_dir["input_fasta_primer"]
    output_file = test_dir["output_primer"]
    tmp_file = test_dir["tmp_file"]
    N = 5  # Number of 'N' characters to use as a separator
    script_path = mock_script

    # Run the merge function
    merge(algo, input_file, output_file, tmp_file, N, script_path)

    # Verify the output file was created
    assert os.path.exists(output_file), "Output file not created"

    # Verify the content of the output file
    with open(output_file, "r") as f:
        output_content = f.read()
    expected_output = ">seq1\nATGCNNNNNCGTA\n>seq2\nTTAANNNNNGGCC\n"
    assert output_content == expected_output, "Output content does not match expected result"

def test_merge_fish_algorithm(test_dir, mock_script):
    """Test the merge function with the 'FISH' algorithm."""
    # Define input parameters
    algo = "FISH"
    input_file = test_dir["input_fasta_fish"]
    output_file = test_dir["output_fish"]
    tmp_file = test_dir["tmp_file"]
    N = 5  # Number of 'N' characters (not used in FISH algorithm)
    script_path = mock_script

    # Run the merge function
    merge(algo, input_file, output_file, tmp_file, N, script_path)

    # Verify the output file was created
    assert os.path.exists(output_file), "Output file not created"

    # Verify the content of the output file
    with open(output_file, "r") as f:
        output_content = f.read()
    expected_output = ">seq1\nATGC\n>seq2\nCGTA\n"
    assert output_content == expected_output, "Output content does not match expected result"

def test_merge_invalid_algorithm(test_dir, mock_script):
    """Test the merge function with an invalid algorithm."""
    # Define input parameters
    algo = "invalid_algorithm"
    input_file = test_dir["input_fasta_primer"]
    output_file = test_dir["output_primer"]
    tmp_file = test_dir["tmp_file"]
    N = 5
    script_path = mock_script

    # Verify that the function raises a KeyError for an invalid algorithm
    with pytest.raises(KeyError, match="Algorithm not implemented yet; see docs or help."):
        merge(algo, input_file, output_file, tmp_file, N, script_path)