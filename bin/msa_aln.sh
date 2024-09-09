#!/bin/bash

# Create a new directory called 'msa'
mkdir -p msa

# Combine all consensus files into one file in the 'msa' directory
cat ../Consensus/*consensus.cov5support80basequal0mask0.fas > msa/combined.fasta

# Display all header lines of the combined fasta file
grep '>' msa/combined.fasta

# Create a multiple sequence alignment using MAFFT
mafft msa/combined.fasta > msa/combined.aln

echo "Script executed successfully. The alignment has been saved in 'msa/combined.aln'."

# INCLUDING PREVIOUSLY PUBLISHED DOGS

cat msa/combined.fasta /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/Canis_latrans.fasta /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/science_dogs_all.with_haps.fasta > msa/combined_with_pub.fasta

# create a msa using MAFFT with previously published dogs
mafft msa/combined_with_pub.fasta > msa/combined_with_pub.aln

echo "Script executed successfully. The alignment has been saved in 'msa/combined_with_pub.aln'."
