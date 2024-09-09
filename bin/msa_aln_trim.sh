#!/bin/bash


# Combine all consensus files into one file in the 'msa' directory
cat ../Consensus/mask/trim/*fas > msa/combined_trim.fasta

# Display all header lines of the combined fasta file
grep '>' msa/combined_trim.fasta

# Create a multiple sequence alignment using MAFFT
mafft msa/combined_trim.fasta > msa/combined_trim.aln

echo "Script executed successfully. The alignment has been saved in 'msa/combined_trim.aln'."

# INCLUDING PREVIOUSLY PUBLISHED DOGS

cat msa/combined_trim.fasta /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/Canis_latrans.fasta /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/science_dogs_all.with_haps.fasta > msa/combined_trim_with_pub.fasta

# create a msa using MAFFT with previously published dogs
mafft msa/combined_trim_with_pub.fasta > msa/combined_trim_with_pub.aln

echo "Script executed successfully. The alignment has been saved in 'msa/combined_trim_with_pub.aln'."
