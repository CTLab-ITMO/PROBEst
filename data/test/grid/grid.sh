# Function to prepare BLAST databases and contigs table
prepare () {
    echo -e "---- Preparation stage ----\n"

    # Clean up previous output directories
    rm -rf data/test/grid/contigs
    rm -rf data/test/grid/blastn_base
    rm -rf data/test/grid/output
    
    # Prepare BLAST database and contigs table for the "true base" dataset
    bash scripts/generator/prep_db.sh \
        -n data/test/grid/blastn_base/true_base \
        -c data/test/grid/contigs \
        -t data/test/grid/.tmp \
        data/test/grid/fasta_base/true_base/*
    
    # Prepare BLAST database and contigs table for the first "false base" dataset
    bash scripts/generator/prep_db.sh \
        -n data/test/grid/blastn_base/false_base \
        -c data/test/grid/contigs \
        -t data/test/grid/tmp \
        data/test/grid/fasta_base/false_base/*
}

# Execute the preparation function
prepare

echo -e "All BLASTn database created" 

# Execute the Python pipeline for probe generation and analysis
python pipeline.py \
    -i data/test/grid/fasta_base/reference.fna \
    -o data/test/grid/output \
    -tb data/test/grid/blastn_base/true_base \
    -fb data/test/grid/blastn_base/false_base \
    -c data/test/grid/contigs \
    -a FISH \
    --PRIMER_PICK_PRIMER 10 \
    --PRIMER_NUM_RETURN 10

# Run grid
python test_parameters.py \
    -p data/test/grid/param_grid_light.json
    -o data/test/grid/out \
    -t 32 -kf