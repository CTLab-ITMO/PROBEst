#!/bin/bash

# PATHS AND EXECUTABLES

#PROBEST=/mnt/c/Users/PC/Desktop/materials/ITMO/thesis/probest/PROBEst/
#BLAST_PARSER=/mnt/c/Users/PC/Desktop/materials/ITMO/thesis/OligoMiner/benchmark/parse_blast_stats.py

BENCH_DATA="$(realpath "$1")"
NUM_ITERATIONS=${2:-10}
TARGETS_DIR=${BENCH_DATA}/target/
OFFTARGETS_DIR=${BENCH_DATA}/offtarget/


# Prepare DB

cat $TARGETS_DIR/*target_* > target_db.fasta
cat $OFFTARGETS_DIR/*off* > offtarget_db.fasta

cat target_db.fasta offtarget_db.fasta > combined.fasta
makeblastdb -in combined.fasta -dbtype nucl -out combined_db

rm ./contigs

$PROBEST/scripts/generator/prep_db.sh \
	-n ./target_db \
	-c ./contigs \
	-t ./db_tmp \
	$TARGETS_DIR/*.fasta

$PROBEST/scripts/generator/prep_db.sh \
        -n ./offtarget_db \
        -c ./contigs \
        -t ./db_tmp \
        $OFFTARGETS_DIR/*.fasta

# Run PROBEST

python $PROBEST/pipeline.py \
	--input $TARGETS_DIR/target_0.fasta \
	--output ./results \
	--true_base ./target_db \
	--false_base ./offtarget_db \
	--contig_table ./contigs \
	--algorithm FISH \
	--iterations $NUM_ITERATIONS \
	--top 50 \
	--mutation_rate 0.05 \
	--set_size 50 \
	--append True

blastn \
        -query ./results/output.fa \
        -db combined_db \
        -out probes_vs_db.blast \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' \
        -task blastn-short

python $BLAST_PARSER ./results/output.fa ./probes_vs_db.blast
