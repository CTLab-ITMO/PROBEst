#installation
#pip install ncbi-genome-download

#file management
kill_structure() {
    for dir in $(ls); do
        if [ -f $dir/*fna.gz ]; then  
            echo $dir  
            fname=$(ls $dir/*fna.gz)
            bn=$(basename $fname)
            mv $fname ./../../$bn
        fi
    done
}

mkdir -p data/test/grid/fasta_base && cd data/test/grid/fasta_base

# Collect true-base
mkdir -p true_base && cd true_base
ncbi-genome-download --formats fasta -P --flat-out -T 562 -s refseq -l complete bacteria 
cd refseq/bacteria
kill_structure
cd ../../
rm -r ./refseq
cd ../

# Collect false-base
mkdir -p false_base && cd false_base
ncbi-genome-download --formats fasta --genera Salmonella -P --flat-out bacteria
cd refseq/bacteria
kill_structure
cd ../../
rm -r ./refseq
cd ../

# Select reference genome for probe generation
cp true_base/GCF_000005845.2_ASM584v2_genomic.fna.gz ./reference.fna.gz
gzip -d ./reference.fna.gz