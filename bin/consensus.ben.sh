#!/bin/bash

# Define directories and reference file
BAM_DIR="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/ChrM/MQ25/dedup"
REF_FILE="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/env_samples/Canis_lupus_familiaris.fasta"
SCRIPT_PATH="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/bin/consensus_from_bam.pl"

# Create the 'Consensus' directory if it doesn't exist
mkdir -p Consensus

# Change to the 'Consensus' directory
cd Consensus

# Loop through all BAM files in the specified directory
for bam_file in "$BAM_DIR"/*.bam; do
    # Extract the base name of the file (without path and extension)
    base=$(basename "$bam_file" .bam)
    
    echo "Processing $base.bam..."

    echo "$SCRIPT_PATH -ref $REF_FILE $bam_file"

    continue

    # Run the Perl script 'consensus_from_bam.pl'
    "$SCRIPT_PATH" \
      -ref "$REF_FILE" \
      "$bam_file"

    echo "Consensus creation process completed for $base.bam"
    echo "----------------------------------------"
done

echo "All BAM files have been processed."
