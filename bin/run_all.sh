##!/bin/bash

#define base path, need to be change to your path
export BASE_PATH="/mnt/expressions/michelle_erkenberg/github/mpi_dogs" 
#export BASE_PATH="/home/michelle/github/mpi_dogs" 

#processing the data, filtering for ChrM and Quality and counting

	#filter for the mitochondrial Chromosom, mapquaility 25% and dedup
#	bash processing/filter_ChrM.sh
#	bash processing/filter_MQ25.sh
#	bash processing/filter_dedup.sh
	
	#count the sequences and for orignial data count also for sequences per chromosom
#	bash processing/count_sequences.sh
#	bash processing/count_sequences_per_chromosom.sh
	
#Call a consensus sequence for each dog using Matthias' perl script
#	bash consensus/consensus.sh	

#msa for the created consensus sequences (renaming the sequence as part of the masking process)

	#copies and renames the consensus sequences for all of our dogs, Undetermined is deleted in the process
#	bash msa/mask.sh

	#cutoff the sequence of all our dogs after TTTTAGG/AAG
#	bash msa/cutoff_seq.sh
	
	#all of your dogs + previously published ones + reference dog combined in one script each + msa
#	bash msa/msa.sh 

#create phylogenetic tree for our dogs and the published dogs with a reference one
#	bash tree/tree.sh

#finding differences and replace the n's (also see diff.sh for more details); it takes the latest aln so far and trims it to 16138bp; afterwards closely related dog sequences are extracted (manually) and differences are highlighted in an csv file; the n's from our dogs are than replaced (if possible) with the bases from the replated dogs; that leads to an csv file with replaced n's and in the long term to an fasta file with out dogs + replaced possions 
#	bash diff/diff.sh

#using the ref dog genome for genome coordinates 
#mkdir -p "$BASE_PATH/data/dog_samples/ref"
#python3 ref/ref_coor.py "$BASE_PATH/data/dog_samples/diff/replaced_seq.related_n.mpi_dogs.added_ref.aln" "$BASE_PATH/data/dog_samples/ref/ref_coordinates.csv" NC_002008.4 fasta

#dogs were living in different offices, this script creates files with just the dogs that lived together 
#bash ref/run_extract_dog.sh 

#finding the private position for our dogs in the environment data
#old version
#python3 env_bam/bam_finder_old.py "$BASE_PATH/data/dog_samples/ref/office_1/5dogs.Heidi.csv" "$BASE_PATH/data/env_samples/Canidae/sample_10.Canidae.Canis_lupus_famili"$BASE_PATH/data/dog_samples/env_bam/env_bam" "Heidi"
#python3 env_bam/bam_finder_new.py "$BASE_PATH/data/dog_samples/ref/office_1/5dogs.Heidi.csv" "$BASE_PATH/data/env_samples/Canidae" "$BASE_PATH/data/dog_samples/env_bam/all_env_Heidi.csv" "Heidi"

bash env_bam/run_bam.sh
