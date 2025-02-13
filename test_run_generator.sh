#create dbs and contigs table
prepare () {
    rm -rf data/test/general/contigs
    rm -rf data/test/general/blastn_base
    rm -rf data/test/general/output
    
    bash scripts/generator/prep_db.sh \
    -n data/test/general/blastn_base/true_base \
    -c data/test/general/contigs \
    -t data/test/general/.tmp \
    data/test/general/fasta_base/true_base/*
    
    bash scripts/generator/prep_db.sh \
    -n data/test/general/blastn_base/false_base_1 \
    -c data/test/general/contigs \
    -t data/test/general/tmp \
    data/test/general/fasta_base/false_base_1/*
    
    bash scripts/generator/prep_db.sh \
    -n data/test/general/blastn_base/false_base_2 \
    -c data/test/general/contigs \
    -t data/test/general/tmp \
    data/test/general/fasta_base/false_base_2/*
}

prepare

#exec
python pipeline.py \
    -i data/test/general/test.fna \
    -o data/test/general/output \
    -tb data/test/general/blastn_base/true_base \
    -fb data/test/general/blastn_base/false_base_1 \
    data/test/general/blastn_base/false_base_2 \
    -c data/test/general/contigs \
    -a FISH \
    --PRIMER_PICK_PRIMER 1 \
    --PRIMER_NUM_RETURN 1 \
    --primer3 primer3_core 

