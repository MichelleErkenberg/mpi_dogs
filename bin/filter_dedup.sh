#!/bin/bash

# Create the dedup subdirectory inside ChrM/MQ25 if it doesn't exist
mkdir -p ChrM/MQ25/dedup

# Loop through all MQ25*.bam files in the ChrM/MQ25 directory
for bam_file in ChrM/MQ25/*.bam
do
    # Extract the filename without path and extension
    filename=$(basename "$bam_file" .bam)
    
    # Deduplicate the BAM file and save in the dedup directory
    samtools rmdup "$bam_file" "ChrM/MQ25/dedup/${filename}_dedup.bam"
    
    # Index the new deduplicated BAM file
    samtools index "ChrM/MQ25/dedup/${filename}_dedup.bam"
    
    echo "Processed and deduplicated: $bam_file"
done

echo "All MQ25 files have been deduplicated and saved in ChrM/MQ25/dedup."
