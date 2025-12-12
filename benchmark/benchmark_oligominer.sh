#!/bin/bash
# Install OligoMiner following instructions on the tool's Github page. (remove NUPACK from environment.yaml)
# Additionally install: 
# - bedtools
# - blast

# PYTHON3 deps
# - biopython
# - pandas

# PATHS AND EXECUTABLES

#PYTHON3=/home/musha/miniconda3/envs/multip/bin/python3
# OLIGOMINER=/mnt/c/Users/PC/Desktop/materials/ITMO/thesis/OligoMiner/
#BLAST_PARSER=/mnt/c/Users/PC/Desktop/materials/ITMO/thesis/OligoMiner/benchmark/parse_blast_stats.py


BENCH_DATA=$(realpath $1)
TARGETS_DIR=${BENCH_DATA}/target/
OFFTARGETS_DIR=${BENCH_DATA}/offtarget/


# PREPARE FILE WITH INSERT COORDINATES

tail -n+2 $BENCH_DATA/metadata.tsv  | awk -F'\t' 'BEGIN {OFS="\t"} {start = $4 - 1; end = start + length($3); print $1, start, end}' | sed 's/-1/0/'  > insert_coords.bed


# MAKE BLAST DB

cat $TARGETS_DIR/*target_* > target_db.fasta
cat $OFFTARGETS_DIR/*off* > offtarget_db.fasta

cat target_db.fasta offtarget_db.fasta > combined.fasta
makeblastdb -in combined.fasta -dbtype nucl -out combined_db

# RUN OligoMiner FOR EACH TARGET SEQUENCE

mkdir probes
mkdir tmp
cd tmp
	
for f in $TARGETS_DIR/*.fasta
do
	cp $f ./
	fasta=$(basename $f)
	prefix=$(echo $fasta | grep -oP 'target_[0-9]+')
	echo Processing $prefix
	python $OLIGOMINER/blockParse.py -l 25 -T 58 -f $fasta &> /dev/null
	bowtie2-build $fasta $prefix.index &> /dev/null
	bowtie2 -x $prefix.index -U $prefix.fastq --no-hd -k 100 --very-sensitive-local -S aligned_${prefix}.sam &> /dev/null
	$OLIGOMINER/outputClean.py -f aligned_${prefix}.sam  -u &> /dev/null
	bedtools intersect -wa -a ./aligned_${prefix}_probes.bed -b ../insert_coords.bed > ${prefix}_probes_filt.bed
	cp *probes_filt.bed ../probes
	rm *
done

cd ../

cat probes/*.bed | awk '{printf(">%s_%s_%s\n%s\n", $1, $2, $3, $4)}' > probes_final.fasta

# BLAST PRODUCED PROBES

blastn \
	-query ./probes_final.fasta \
	-db combined_db \
	-out probes_vs_db.blast \
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' \
	-task blastn-short


$PYTHON3 $BLAST_PARSER ./probes_final.fasta ./probes_vs_db.blast
