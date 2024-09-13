#!/bin/bash

# Set working directory
WORKDIR="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/Consensus/msa"

# Set output directory
OUTDIR="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/Consensus/msa"

# Change to working directory
cd "$WORKDIR"

# Create phylogenetic tree with MPI dogs (org)
FastTree -nt < combined.aln > "$OUTDIR/mpi_dogs_tree.nwk"

echo "Phylogenetic tree for MPI dogs created and saved as '$OUTDIR/mpi_dogs_tree.nwk'."

# Create phylogenetic tree with MPI dogs and published dogs
FastTree -nt < combined_with_pub.aln > "$OUTDIR/mpi_dogs_with_pub_tree.nwk"

echo "Phylogenetic tree for MPI dogs and published dogs created and saved as '$OUTDIR/mpi_dogs_with_pub_tree.nwk'."

# Create phylogenetic tree with MPI dogs (trim)
FastTree -nt < combined_trim.aln > "$OUTDIR/mpi_dogs_tree_trim.nwk"

echo "Phylogenetic tree for MPI dogs created and saved as '$OUTDIR/mpi_dogs_tree_trim.nwk'."

# Create phylogenetic tree with MPI dogs and published dogs
FastTree -nt < combined_trim_with_pub.aln > "$OUTDIR/mpi_dogs_with_pub_tree_trim.nwk"

echo "Phylogenetic tree for MPI dogs and published dogs created and saved as '$OUTDIR/mpi_dogs_with_pub_tree_trim.nwk'."

# Create phylogenetic tree with MPI dogs (cutoff)
FastTree -nt < combined_cutoff.aln > "$OUTDIR/mpi_dogs_tree_cutoff.nwk"

echo "Phylogenetic tree for MPI dogs created and saved as '$OUTDIR/mpi_dogs_tree_cutoff.nwk'."

# Create phylogenetic tree with MPI dogs and published dogs (cutoff)
FastTree -nt < combined_cutoff_with_pub.aln > "$OUTDIR/mpi_dogs_with_pub_tree_cutoff.nwk"

echo "Phylogenetic tree for MPI dogs and published dogs created and saved as '$OUTDIR/mpi_dogs_with_pub_tree_cutoff.nwk'."

# Create phylogenetic tree with MPI dogs and published dogs (cutoff) and reference dog NC_002008.4
FastTree -nt < combined_cutoff_with_pub.NC_002008.4.aln > "$OUTDIR/mpi_dogs_with_pub_tree_cutoff.NC_002008.4.nwk"

echo "Phylogenetic tree for MPI dogs and published dogs created and saved as '$OUTDIR/mpi_dogs_with_pub_tree_cutoff.NC_002008.4.nwk'."
