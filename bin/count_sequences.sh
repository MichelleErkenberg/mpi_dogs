#!/bin/bash

# Filename: count_sequences.sh

# Use the provided directory or the current one if none is specified
dir="${1:-.}"

# Define the output CSV file
output_file="$dir/sequence_counts.csv"

# Write the header to the CSV file
echo "Filename,SequenceCount" > "$output_file"

# Iterate over all BAM files in the current directory
for bam_file in *.bam
do
  # Count the number of sequences in the BAM file
  sequence_count=$(samtools view -c "$bam_file")
  
  # Print the filename and the number of sequences to the console
  echo "File: $bam_file - Number of sequences: $sequence_count"
  
  # Append the data to the CSV file
  echo "$bam_file,$sequence_count" >> "$output_file"
done

echo "Results have been written to $output_file"
