#!/bin/bash

# Combine all consensus files into one file
cat ../Consensus/*consensus.cov5support80basequal0mask0.fas > combined.fasta

# Display all header lines of the combined fasta file
grep '>' combined.fasta

# Create a multiple sequence alignment using MAFFT
mafft combined.fasta > combined.aln

echo "Script executed successfully. The alignment has been saved in 'combined.aln'."
