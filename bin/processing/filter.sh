#!/bin/bash

# Set base directory 
BASE_PATH="/path/to/your/directory"

# Change to the working directory
cd "$BASE_PATH/data/dog_samples" || exit 1

### Step 1: ChrM Extraction ###
mkdir -p ChrM
echo "Starting ChrM extraction..."
for bam_file in *.bam; do
  base_name=$(basename "$bam_file" .bam)
  samtools view -b "$bam_file" "chrM" > "ChrM/${base_name}_ChrM.bam"
  echo "Processed: $bam_file → ChrM/${base_name}_ChrM.bam"
done
echo "ChrM extraction completed."

### Step 2: MAPQ25 Filtering ###
mkdir -p ChrM/MQ25
echo -e "\nStarting MAPQ25 filtering..."
for bam_file in ChrM/*.bam; do
  filename=$(basename "$bam_file" .bam)
  samtools view -bq 25 "$bam_file" > "ChrM/MQ25/${filename}_MQ25.bam"
  samtools index "ChrM/MQ25/${filename}_MQ25.bam"
  echo "Processed: $bam_file → MQ25/${filename}_MQ25.bam"
done
echo "MAPQ25 filtering completed."

### Step 3: Deduplication ###
mkdir -p ChrM/MQ25/dedup
echo -e "\nStarting deduplication..."
for bam_file in ChrM/MQ25/*.bam; do
  filename=$(basename "$bam_file" .bam)
  samtools rmdup "$bam_file" "ChrM/MQ25/dedup/${filename}_dedup.bam"
  samtools index "ChrM/MQ25/dedup/${filename}_dedup.bam"
  echo "Deduplicated: $bam_file → dedup/${filename}_dedup.bam"
done
echo -e "\nAll processing steps completed successfully!"
