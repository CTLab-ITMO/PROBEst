#create dbs and contigs table
prepare () {
    rm data/test/contigs
    rm -r data/test/blastn_base
    rm -r data/test/output
    
    bash scripts/generator/prep_db.sh \
    -n data/test/blastn_base/true_base \
    -c data/test/contigs \
    -t data/test/.tmp \
    data/test/fasta_base/true_base/*
    
    bash scripts/generator/prep_db.sh \
    -n data/test/blastn_base/false_base_1 \
    -c data/test/contigs \
    -t data/test/tmp \
    data/test/fasta_base/false_base_1/*
    
    bash scripts/generator/prep_db.sh \
    -n data/test/blastn_base/false_base_2 \
    -c data/test/contigs \
    -t data/test/tmp \
    data/test/fasta_base/false_base_2/*
}

prepare

#exec
python pipeline.py \
    -i data/test/test.fna \
    -o data/test/output \
    -tb data/test/blastn_base/true_base \
    -fb data/test/blastn_base/false_base_1 \
    data/test/blastn_base/false_base_2 \
    -c data/test/contigs \
    -a FISH \
    --PRIMER_PICK_PRIMER 1 \
    --PRIMER_NUM_RETURN 1 \
    --primer3 primer3_core 

