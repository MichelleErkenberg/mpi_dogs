#!/bin/bash

# Create phylogenetic tree with MPI dogs (org)
FastTree -nt < combined.aln > mpi_dogs_tree.nwk

echo "Phylogenetic tree for MPI dogs created and saved as 'msa/mpi_dogs_tree.nwk'."

# Create phylogenetic tree with MPI dogs and published dogs
FastTree -nt < combined_with_pub.aln > mpi_dogs_with_pub_tree.nwk

echo "Phylogenetic tree for MPI dogs and published dogs created and saved as 'msa/mpi_dogs_with_pub_tree.nwk'."

# Create phylogenetic tree with MPI dogs (trim)
FastTree -nt < combined_trim.aln > mpi_dogs_tree_trim.nwk

echo "Phylogenetic tree for MPI dogs created and saved as 'msa/mpi_dogs_tree_trim.nwk'."

# Create phylogenetic tree with MPI dogs and published dogs
FastTree -nt < combined_trim_with_pub.aln > mpi_dogs_with_pub_tree_trim.nwk

echo "Phylogenetic tree for MPI dogs and published dogs created and saved as 'msa/mpi_dogs_with_pub_tree_trim.nwk'."

# Create phylogenetic tree with MPI dogs (cutoff)
FastTree -nt < combined_cutoff.aln > mpi_dogs_tree_cutoff.nwk

echo "Phylogenetic tree for MPI dogs created and saved as 'msa/mpi_dogs_tree_cutoff.nwk'."

# Create phylogenetic tree with MPI dogs and published dogs (cutoff)
FastTree -nt < combined_cutoff_with_pub.aln > mpi_dogs_with_pub_tree_cutoff.nwk

echo "Phylogenetic tree for MPI dogs and published dogs created and saved as 'msa/mpi_dogs_with_pub_tree_cutoff.nwk'."
