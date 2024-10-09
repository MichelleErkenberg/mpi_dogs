##!/bin/bash

#define base path, need to be change to your path
export BASE_PATH="/mnt/expressions/michelle_erkenberg/github/mpi_dogs" 

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

#highlighting the differences between the dogs and replacing the n's (see diff.sh for details)
	bash diff/diff.sh
