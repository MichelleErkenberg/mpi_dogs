#!/bin/bash
	
#data processing and filtering, getting an overview of the amount and quality of the data	
	#renaming our *.bam files in the mpi dogs
	bash filter_data/rename.sh
	
	#filter the data for mitochondrial Chromosom (ChrM), Mapquality 25% and deduped
	bash filter_data/filter_ChrM.sh 
	bash filter_data/filter_MQ25.sh
	bash filter_data/filter_dedup.sh
	
	#count sequences and sequences per chromosom (just for the original data)
	bash filter_data/count_sequences.sh
	bash filter_data/count_sequences_per_chromosom.sh
	
	
	
