##!/bin/bash

#processing the data: renaming, filtering for ChrM and Quality and counting

	#renaming the *.bam files (in mpi dogs)
	bash processing/rename.sh
	
	#filter for the mitochondrial Chromosom, mapquaility 25% and dedup
	bash processing/filter_ChrM.sh
	bash processing/filter_MQ25.sh
	bash processing/filter_dedup.sh
	
	#count the sequences and for orignial data count also for sequences per chromosom
	bash processing/count_sequences.sh
	bash processing/count_sequences_per_chromosom.sh
	
#Call a consensus sequence for each dog using Matthias' perl script

	bash consensus/consensus.sh	
