#!/bin/bash

#  Combine all consensus files into one file in the 'msa' directory
cat ../Consensus/mask/cutoff/*fas > msa/combined_cutoff.fasta

# Display all header lines of the combined fasta file
grep '>' msa/combined_cutoff.fasta

# Create a multiple sequence alignment using MAFFT
mafft msa/combined_cutoff.fasta > msa/combined_cutoff.aln

echo "Script executed successfully. The alignment has been saved in 'msa/combined_cutoff.aln'."

# INCLUDING PREVIOUSLY PUBLISHED DOGS

cat msa/combined_cutoff.fasta /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/Canis_latrans.fasta /mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/science_dogs/science_dogs_all.with_haps.fasta > msa/combined_cutoff_with_pub.fasta

# create a msa using MAFFT with previously published dogs
mafft msa/combined_cutoff_with_pub.fasta > msa/combined_cutoff_with_pub.aln

echo "Script executed successfully. The alignment has been saved in 'msa/combined_cutoff_with_pub.aln'."
