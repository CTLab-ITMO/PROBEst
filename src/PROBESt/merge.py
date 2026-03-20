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


import subprocess
from shutil import copyfile


def _primer_fasta_is_paired(fasta_path: str) -> bool:
    """True if Primer3-style LEFT/RIGHT records; False if already merged or mutation output."""
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                return "_LEFT" in line or "_RIGHT" in line
    return False


def _parse_primer_side_header(line: str, side: str) -> str:
    """Id for one amplicon: full header without *_LEFT / *_RIGHT (keeps primer pair index)."""
    s = line.rstrip("\r\n")
    suf = "_LEFT" if side == "LEFT" else "_RIGHT"
    if not s.endswith(suf):
        raise LookupError(f"Malformed {side} primer header: {s!r}")
    return s[: -len(suf)]


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

    if algo == "primer" and not _primer_fasta_is_paired(input):
        copyfile(input, output)
        print("Primers merged successfully.")
        return

    fasta2table_command = f"bash {script_path}fasta2table.sh -i {input} -o {tmp}"
    subprocess.run(fasta2table_command, shell=True)

    if algo == "primer":
        with open(tmp, "r") as fasta_table, open(output, "w") as output_fasta:
            sepN = "N" * NNN

            for i, line in enumerate(fasta_table):
                dev = i % 4
                if dev == 0:
                    inp_name = _parse_primer_side_header(line, "LEFT")
                elif dev == 1:
                    left = line.strip()
                elif dev == 2:
                    right_base = _parse_primer_side_header(line, "RIGHT")
                    if inp_name != right_base:
                        raise LookupError(f"Error on line {i}: check names in FASTA file")
                else:
                    right = line.strip()
                    output_fasta.write(f"{inp_name}\n{left}{sepN}{right}\n")

        print("Primers merged successfully.")

    elif algo == "FISH":
        copyfile(tmp, output)

    else:
        raise KeyError("Algorithm not implemented yet; see docs or help.")
