#!/bin/bash

# Combine all consensus files into one file
cat ../Consensus/*consensus.cov5support80basequal0mask0.fas > combined.fasta

# Display all header lines of the combined fasta file
grep '>' combined.fasta

# Create a multiple sequence alignment using MAFFT
mafft combined.fasta > combined.aln

echo "Script executed successfully. The alignment has been saved in 'combined.aln'."

# INCLUDING PREVIOUSLY PUBLISHED DOGS

cat combined.fasta /mnt/scratch/praktikum2024_mol_anth/data/Canis_latrans.fasta /mnt/scratch/praktikum2024_mol_anth/data/science_dogs/science_dogs_all.with_haps.fasta > combined_with_pub.fasta

# create a msa using MAFFT with previously published dogs
mafft combined_with_pub.fasta > combined_with_pub.aln

echo "Script executed successfully. The alignment has been saved in 'combined_with_pub.aln'."
