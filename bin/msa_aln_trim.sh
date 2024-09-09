#!/bin/bash

# Set the working directory
work_dir="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/Consensus/mask/trim"
msa_dir="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/Consensus/msa"

# Change to the working directory
cd "$work_dir" || exit 1

# Ensure the MSA directory exists
mkdir -p "$msa_dir"

# Combine all consensus files into one file in the 'msa' directory
cat *.fas > "$msa_dir/combined_trim.fasta"

# Display all header lines of the combined fasta file
grep '>' "$msa_dir/combined_trim.fasta"

# Create a multiple sequence alignment using MAFFT
mafft "$msa_dir/combined_trim.fasta" > "$msa_dir/combined_trim.aln"

echo "First part executed successfully. The alignment has been saved in '$msa_dir/combined_trim.aln'."

# INCLUDING PREVIOUSLY PUBLISHED DOGS

cat "$msa_dir/combined_trim.fasta"  /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/Canis_latrans.fasta /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/science_dogs_all.with_haps.fasta > "$msa_dir/combined_trim_with_pub.fasta"

# create a msa using MAFFT with previously published dogs
mafft "$msa_dir/combined_trim_with_pub.fasta" > "$msa_dir/combined_trim_with_pub.aln"

echo "Script executed successfully. The alignment has been saved in '$msa_dir/combined_trim_with_pub.aln'."
