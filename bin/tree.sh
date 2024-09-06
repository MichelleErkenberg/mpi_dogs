#!/bin/bash

# Create phylogenetic tree with MPI dogs
FastTree -nt < msa/combined.aln > msa/mpi_dogs_tree.nwk

echo "Phylogenetic tree for MPI dogs created and saved as 'msa/mpi_dogs_tree.nwk'."

# Create phylogenetic tree with MPI dogs and published dogs
FastTree -nt < msa/combined_with_pub.aln > msa/mpi_dogs_with_pub_tree.nwk

echo "Phylogenetic tree for MPI dogs and published dogs created and saved as 'msa/mpi_dogs_with_pub_tree.nwk'."
