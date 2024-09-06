#!/bin/bash

# Create output directory
mkdir -p ChrM

# Iterate over all BAM files in the current directory
for bam_file in *.bam; do
  # Get the filename without extension
  base_name=$(basename "$bam_file" .bam)
  
  # Output file in the chrM directory
  output_file="ChrM/${base_name}_ChrM.bam"
  
  # Filter chrM and write to the new file
  samtools view -b "$bam_file" "chrM" > "$output_file"
  
  echo "Processed $bam_file -> $output_file"
done

echo "All files have been processed and saved in the ChrM directory."
