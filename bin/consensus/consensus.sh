#!/bin/bash

# Define directories and reference file
BAM_DIR="../data/dog_samples/ChrM/MQ25/dedup"
REF_FILE="../data/science_dogs/Canis_lupus_familiaris.fasta"
SCRIPT_PATH="/consensus/consensus_from_bam.pl"

# Create the 'consensus' directory if it doesn't exist
mkdir -p ../data/dog_samples/consensus

# Change to the 'consensus' directory
cd ../data/dog_samples/consensus

# Loop through all BAM files in the specified directory
for bam_file in "$BAM_DIR"/*.bam; do
    # Extract the base name of the file (without path and extension)
    base=$(basename "$bam_file" .bam)
    
    echo "Processing $base.bam..."

    # Run the Perl script 'consensus_from_bam.pl'
    "$SCRIPT_PATH" \
      -ref "$REF_FILE" \
      "$bam_file"

    echo "Consensus creation process completed for $base.bam"
    echo "----------------------------------------"
done

echo "All BAM files have been processed."
