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

def uniline_fasta(args, out):    
    uniline = "bash " + args.script_path + "uniline_fa.sh"
    uniline += " -i " + args.input
    uniline += " -o " + out + "input.fa"
    subprocess.run(uniline, shell=True)

def blastn_function(args: list):
    blastn_comand = args.blastn + " -num_threads " + \
        args.threads + " -outfmt '6 qseqid sseqid evalue sstart send ppos mismatch' " + \
        " -word_size " + args.word_size + \
        " -reward " + args.reward + \
        " -penalty " + args.penalty + \
        " -gapopen " + args.gapopen + \
        " -gapextend " + args.gapextend + \
        " -evalue " + args.evalue

    return blastn_comand


def probe_check_function(args: list):
    probe_check_comand = "bash " + args.script_path + "/probe_check.sh" + \
        " -p " + args.script_path + "/probe_filt.py" + \
        " -d " + args.contig_table + \
        " -m " + str(args.top) + \
        " -e " + args.max_mismatch + \
        " -i " + args.min_ident + \
        " -a " + args.multimap_max + \
        " -b " + args.negative_max
    
    return probe_check_comand
