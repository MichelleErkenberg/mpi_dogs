##!/bin/bash

#change to directory
bam_files="$BASE_PATH/data/dog_samples/processing/"

# Function to process BAM files
process_bam_file() {
  local bam_file="$1"
  local current_dir="$(dirname "$bam_file")"
  local output_file="$current_dir/$(basename "$bam_file" .bam)_sequence_counts.csv"

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
}

# Recursively iterate through all directories and process BAM files
find "$bam_files" -type f -name "*.bam" | while read -r bam_file; do
  process_bam_file "$bam_file"
done
