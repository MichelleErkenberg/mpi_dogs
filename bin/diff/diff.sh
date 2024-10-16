##!/bin/bash

TREE_DIR="$BASE_PATH/data/dog_samples/tree"
OUTDIR="$BASE_PATH/data/dog_samples/diff"

#create OUTDIR if it doesn't exist
mkdir -p "$OUTDIR"

#create and trim the .aln file to cut all sequences to a defined length of 16138 bp
bash diff/trim_msa.sh

#Extract the closely related reference dogs for our mpi dogs (if possible) + alignment; mpi_dogs replaced by mpi_dogs are not used for comparison
python3 diff/extract_sequences.py "$TREE_DIR/mpi_dogs_pub.ref.nwk" "$OUTDIR/combined_pub.ref.trimmed.aln" names.txt "$OUTDIR/related_seq.fasta"
mafft "$OUTDIR/related_seq.fasta" > "$OUTDIR/related_seq.aln"

#creating .csv-files to highlight variations between the dogs
python3 diff/dog_vari.py "$OUTDIR/related_seq.aln" "$OUTDIR/sequence_diff.csv"

#creates a .csv-file where N's are replaced by bases from closely related dogs (if they are the same)
python3 diff/replace_n.py names.txt "$OUTDIR/sequence_diff.csv" "$OUTDIR/replaced_n.csv"

#uses the replaced_n.csv to create a new .aln-file with replaced n's
python3 diff/replaced_aln.py "$OUTDIR/replaced_n.csv" "$OUTDIR/related_seq.aln" "$OUTDIR/replaced_n.aln"  

#write the mpi dog sequences with replaced n's into a new fasta 
bash diff/mpi_fasta.sh
