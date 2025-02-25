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
