#!/bin/bash

# PATHS AND EXECUTABLES

#BENCH_DATA=/mnt/c/Users/PC/Desktop/materials/ITMO/thesis/OligoMiner/benchmark/test_data/
#BLAST_PARSER=/mnt/c/Users/PC/Desktop/materials/ITMO/thesis/OligoMiner/benchmark/parse_blast_stats.py

BENCH_DATA=$(realpath $1)
TARGETS_DIR=${BENCH_DATA}/target/
OFFTARGETS_DIR=${BENCH_DATA}/offtarget/


tail -n+2 $BENCH_DATA/metadata.tsv  | awk -F'\t' 'BEGIN {OFS="\t"} {start = $4 - 1; end = start + length($3); print $1, start, end}' | sed 's/-1/0/'  > insert_coords.bed

# MAKE BLAST DB

cat $TARGETS_DIR/*target_* > target_db.fasta
cat $OFFTARGETS_DIR/*off* > offtarget_db.fasta

cat target_db.fasta offtarget_db.fasta > combined.fasta
makeblastdb -in combined.fasta -dbtype nucl -out combined_db

# RUN primer3

mkdir probes
mkdir tmp
cd tmp

for f in $TARGETS_DIR/*.fasta
do
    cp $f ./
    fasta=$(basename $f)
    prefix=$(echo $fasta | grep -oP 'target_[0-9]+')
    python $PRIMER3_PARSER $fasta
    primer3_core ./primer3_template --output primer3.output
    cat primer3.output \
        | grep -P 'PRIMER_RIGHT_[0-9]+_SEQUENCE' \
        | sed "s/^/>${prefix}_/" \
        | tr '=' '\n' > ../probes/${prefix}_probes.fasta
    rm *
done

cd ../
cat probes/*.fasta > probes_final.fasta

# BLAST PRODUCED PROBES

blastn \
	-query ./probes_final.fasta \
	-db combined_db \
	-out probes_vs_db.blast \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' \
	-task blastn-short


python $BLAST_PARSER ./probes_final.fasta ./probes_vs_db.blast