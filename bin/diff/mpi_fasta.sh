#!/bin/bash

# Define environment
DIR="$BASE_PATH/data/dog_samples/diff"
aln_file="$DIR/replaced_n.aln"
names_file="$BASE_PATH/bin/names.txt"

# Check if files exist
if [ ! -f "$aln_file" ] || [ ! -f "$names_file" ]; then
    echo "Error: $aln_file or $names_file does not exist."
    exit 1
fi

# Create a temporary file for output
temp_output=$(mktemp)

# Ensure temporary file is removed on script exit
trap 'rm -f "$temp_output"' EXIT

# Read names.txt into array
mapfile -t names < "$names_file"

# Read alignment file and find dog names
extract_sequence=false
while IFS= read -r line; do
    if [[ $line =~ ^'>' ]]; then
        sequence_name=${line#>}
        if [[ " ${names[@]} " =~ " ${sequence_name} " ]]; then
            extract_sequence=true
            echo "$line" >> "$temp_output"
        else
            extract_sequence=false
        fi
    elif [ "$extract_sequence" = true ]; then
        echo "$line" >> "$temp_output"
    fi
done < "$aln_file"

# Save results to output file
output_file="$DIR/mpi_dogs_replaced_n.fasta"
mv "$temp_output" "$output_file"

echo "MPI dog sequences saved to $output_file"
