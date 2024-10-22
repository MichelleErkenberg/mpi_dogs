#!/bin/bash

# Define environment
DIR="$BASE_PATH/data/dog_samples/diff"
aln_file="$DIR/replaced_n.aln"
names_file="$BASE_PATH/bin/names.txt"
output_file="$DIR/mpi_dogs_replaced_n.fasta"

# Create a grep pattern from names.txt
grep_pattern=$(sed 's/^/^>/; s/$/|/' "$names_file" | tr -d '\n' | sed 's/|$//')

# Extract sequences and save to output file
grep -A 1 -f <(echo "$grep_pattern") "$aln_file" | sed '/^--$/d' > "$output_file"

echo "MPI dog sequences saved to $output_file"

