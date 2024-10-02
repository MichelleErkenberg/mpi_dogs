#!/bin/bash

#change to directory
cd ../data/dog_samples

# Iterate over all BAM files in the specified directory
for bam_file in *.bam; do
  # Define the output file for each BAM file
  output_file="$(basename "$bam_file" .bam)_sequence_counts.csv"

  # Write the header to the CSV file
  echo "Chromosome,SequenceCount" > "$output_file"

  # Use samtools idxstats to get the count of sequences per chromosome
  samtools idxstats "$bam_file" | while read -r line
  do
    # Extract chromosome name and sequence count
    chromosome=$(echo "$line" | cut -f 1)
    sequence_count=$(echo "$line" | cut -f 3)

    # Print the filename, chromosome, and the number of sequences
    echo "File: $(basename "$bam_file") - Chromosome: $chromosome - Number of sequences: $sequence_count"

    # Append the data to the CSV file
    echo "$chromosome,$sequence_count" >> "$output_file"
  done

  echo "Results for $(basename "$bam_file") have been written to $output_file"
done
