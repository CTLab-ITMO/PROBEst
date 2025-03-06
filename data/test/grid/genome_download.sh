#installation
#pip install ncbi-genome-download

#file management
kill_structure() {
    for dir in $(ls); do
        fname=$(ls $dir/*fna.gz)
        mv $fname ./../../
    done
}

mkdir -p data/test/grid/fasta_base && cd data/test/grid/fasta_base

# Collect true-base
mkdir -p true_base && cd true_base
ncbi-genome-download --formats fasta --genera Borrelia bacteria
cd refseq/bacteria
kill_structure
cd ../../
rm -r ./refseq
cd ../

# Collect false-base
mkdir -p false_base && cd false_base
ncbi-genome-download --formats fasta --genera Rickettsia bacteria
cd refseq/bacteria
kill_structure
cd ../../
rm -r ./refseq
cd ../

# Select reference genome for probe generation
cp true_base/GCF_000012065.2_ASM1206v2_genomic.fna.gz ./reference.fna.gz
gzip -d ./reference.fna.gz