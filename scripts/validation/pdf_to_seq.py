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


from pdfminer.high_level import extract_text
import re
import argparse
import sys
import os
import pathlib

def find_probes(path_to_file, min_length, atcg_threshold):
    # Extract text from the PDF
    text = extract_text(path_to_file)

    # Regex pattern to match DNA probe candidates
    pattern = r"(?:5'?\s*)?(?:[A-Za-z0-9]+-)?[ATCG]+(?:-[A-Za-z0-9]+)?(?:\s*3'?)?"
    if atcg_threshold == 1:
        pattern = r"[^ATCG]*([ATCG]+)[^ATCG]*"

    matches = re.findall(pattern, text, re.VERBOSE | re.IGNORECASE)

    probes = []
    for candidate in matches:
        # Clean candidate (remove extra spaces)
        candidate = re.sub(r"\s+", "", candidate.strip())
        
        # Count A/T/C/G characters (case-insensitive)
        atcg_count = sum(1 for c in candidate.upper() if c in {"A", "T", "C", "G"})
        total_chars = len(candidate)
        
        # Skip if too short or low ATCG content
        if total_chars == 0 or atcg_count < min_length:
            continue
        if (atcg_count / total_chars) < atcg_threshold:
            continue
        
        probes.append(candidate)

    return probes


def main(cmdline=None):
    if cmdline is None:
        cmdline=sys.argv[1:]
    parser = argparse.ArgumentParser(
        description='This scripts validate JSON from LLM against determined schema with data types and RegExp')
    parser.add_argument('file_path', type=str, help=
                        "Path to article PDF file")
    parser.add_argument('--min_len','-m', default=9, type=int, help=
                    "Minimal length of probe sequence, defaults to 9")
    parser.add_argument('--atcg_threshold','-t', default=0.6, type=float, help=
                "Minimal portion of characters ATCG in probe sequence (probe can contain modifications or 5' and 3' symbols),\
                setting this parameter to 1 will make script extract only ATCG containing sequences, defaults to 0.6")
    parser.add_argument('--output','-o', default='probes.txt', type=str, help=
            "Name of output file, default is probes.txt")
    
    args = parser.parse_args(cmdline)
    print(args)
    file_path = args.file_path
    if not os.path.isabs(file_path):
            file_path = pathlib.Path(os.getcwd()) / pathlib.Path(file_path)

    result = find_probes(file_path, int(args.min_len), float(args.atcg_threshold))
    if len(result) > 0:
        print("DNA Probes Found:", result)
        with open(args.output, 'w') as f:
            for p in result:
                f.write(p)
                f.write('\n')
    else:
        print('No probes found, check PDF file')

if __name__ == '__main__':
    main()