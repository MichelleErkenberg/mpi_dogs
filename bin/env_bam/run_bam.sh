#!/bin/bash

#creating and define new env_bam file as csv output file
mkdir -p "$BASE_PATH/data/dog_samples/env_bam"
env_bam="$BASE_PATH/data/dog_samples/env_bam"

#creating *.bam.bai data for every bam to process
for bam in "$BASE_PATH/data/env_samples/Canidae/*.bam"
do
	if [ ! -f "${bam}.bai" ]; then
        echo "Indexing $bam..."
        samtools index "$bam"
	else
        echo "Index for $bam already exists, skipping."
    	fi
done

echo "BAM indexing process completed. All BAM files are now indexed or were already indexed."

#creating seperated file for our dogs to structure
mkdir -p "$env_bam/Heidi"

#finding 
python3 bam_finder.py "$BASE_PATH/data/dog_samples/ref/office_1/5dogs.Heidi.csv" "$BASE_PATH/data/env_samples/Canidae" "$env_bam/Heidi" "Heidi"
