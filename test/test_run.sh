##!/bin/bash

#align the fasta file
mafft data/test.fasta > data/test.aln

#write the aln file into a csv file to replace in the n's afterwards in this file
python3 aln_to_csv.py data/test.aln data/test.csv

#replacing the n's from our dogs with the references
python3 replace_n.py names.txt data/test.csv data/n_re.csv

#write the csv file back again into a fasta format
python3 replaced_n_fasta.py data/n_re.csv data/replaced_n.fasta

#write the mpi dog sequences with replaced n's into a new fasta 
seqkit grep -f names.txt data/replaced_n.fasta > data/mpi_replaced_n.fasta

#add the ref dog file
cat data/ref.fasta data/mpi_replaced_n.fasta > data/mpi_ref.fasta

#align ref and mpi dog
mafft data/mpi_ref.fasta > data/mpi_ref.aln

#references dog coordinates
python3 ref_coor.py data/mpi_ref.aln data/ref_coordinates.csv Reference fasta




