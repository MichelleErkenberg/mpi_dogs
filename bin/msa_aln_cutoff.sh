#!/bin/bash

# Set the working directory
work_dir="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/Consensus/mask/cutoff"
msa_dir="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/Consensus/msa"

# Change to the working directory
cd "$work_dir" || exit 1

# Ensure the MSA directory exists
mkdir -p "$msa_dir"

#  Combine all consensus files into one file in the 'msa' directory
cat *.fas > "$msa_dir/combined_cutoff.fasta"

# Display all header lines of the combined fasta file
grep '>' "$msa_dir/combined_cutoff.fasta"

# Create a multiple sequence alignment using MAFFT
mafft "$msa_dir/combined_cutoff.fasta" > "$msa_dir/combined_cutoff.aln"

echo "Script executed successfully. The alignment has been saved in 'msa/combined_cutoff.aln'."

# INCLUDING PREVIOUSLY PUBLISHED DOGS

cat "$msa_dir/combined_cutoff.fasta" /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/Canis_latrans.fasta /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/science_dogs_all.with_haps.fasta > "$msa_dir/combined_cutoff_with_pub.fasta"

# create a msa using MAFFT with previously published dogs
mafft "$msa_dir/combined_cutoff_with_pub.fasta" > "$msa_dir/combined_cutoff_with_pub.aln"

echo "Script executed successfully. The alignment has been saved in 'msa/combined_cutoff_with_pub.aln'."

# Including reference dog NC_002008.4

cat "$msa_dir/combined_cutoff_with_pub.aln" /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/Canis_lupus_familiaris.fasta > "$msa_dir/combined_cutoff_with_pub.NC_002008.4.fasta"

#create a msa using MAFFT with reference dog NC_002008.4

mafft "$msa_dir/combined_cutoff_with_pub.NC_002008.4.fasta" > "$msa_dir/combined_cutoff_with_pub.NC_002008.4.aln"

echo "Script executed successfully. The alignment has been saved in 'msa/combined_cutoff_with_pub.NC_002008.4.aln'."

