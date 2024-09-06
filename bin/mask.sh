#!/bin/bash

# File containing the names
NAMES_FILE="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/bin/names.txt"

# FASTA file
FASTA_FILE="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/Consensus/combined_with_pub.fasta"

# Output directory
OUTPUT_DIR="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/Consensus/mask"

# Check if the files exist
if [ ! -f "$NAMES_FILE" ]; then
    echo "Error: The file $NAMES_FILE does not exist."
    exit 1
fi

if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: The file $FASTA_FILE does not exist."
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Read each name from the file and execute seqkit
while IFS= read -r name; do
    echo "Processing: $name"
    seqkit grep -r -p "$name" "$FASTA_FILE" > "${OUTPUT_DIR}/${name}.fas"
done < "$NAMES_FILE"

echo "Done! Output files are in the '$OUTPUT_DIR' directory."
