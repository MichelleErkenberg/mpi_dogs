#!/bin/bash

# Create the MQ25 subdirectory inside ChrM if it doesn't exist
mkdir -p ChrM/MQ25

# Loop through all BAM files in the ChrM directory
for bam_file in ChrM/*.bam
do
    # Extract the filename without path and extension
    filename=$(basename "$bam_file" .bam)
    
    # Filter reads with MAPQ ≥ 25 and save in ChrM/MQ25 directory
    samtools view -bq 25 "$bam_file" > "ChrM/MQ25/${filename}_MQ25.bam"
    
    # Index the new BAM file
    samtools index "ChrM/MQ25/${filename}_MQ25.bam"
    
    echo "Processed: $bam_file"
done

echo "All files have been processed and saved in ChrM/MQ25."
