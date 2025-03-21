#!/bin/bash

# Define the base directory
base_dir="${BASE_PATH}/data/dog_samples/processing"

# Function to extract the sample name from the filename
extract_sample_name() {
    local filename="$1"
    echo "$filename" | sed -n 's/.*s_all_\(.*\)_S.*/\1/p'
}

# Function to process BAM files in the bam_files directory
process_bam_files() {
    local dir="${base_dir}/bam_files"
    local output_file="${base_dir}/bam_files_sequence_counts.csv"

    # Check if there are any BAM files in this directory
    if ! compgen -G "$dir"/*.bam > /dev/null; then
        echo "No BAM files found in $dir"
        return
    fi

    # Initialize the CSV file with headers
    echo "Sample,SequenceCount" > "$output_file"

    # Process each BAM file
    for bam_file in "$dir"/*.bam; do
        sample_name=$(extract_sample_name "$(basename "$bam_file")")
        sequence_count=$(samtools view -c "$bam_file")
        echo "${sample_name},${sequence_count}" >> "$output_file"
        echo "Processed: $sample_name - Count: $sequence_count"
    done

    echo "Results for bam_files have been written to $output_file"
}

# Function to process ChrM and subdirectories
process_chrm_directories() {
    local output_file="${base_dir}/chrm_sequence_counts.csv"

    # Initialize the CSV file with headers
    echo "Sample,ChrM,MQ25,dedup" > "$output_file"

    # Find all sample directories in ChrM
    find "${base_dir}/ChrM" -mindepth 1 -maxdepth 1 -type d | while read -r sample_dir; do
        sample_name=$(extract_sample_name "$(basename "$sample_dir")")
        chrm_count=$(samtools view -c "${sample_dir}/s_all_${sample_name}_S_ChrM.bam" 2>/dev/null || echo "N/A")
        mq25_count=$(samtools view -c "${sample_dir}/MQ25/s_all_${sample_name}_S_ChrM_MQ25.bam" 2>/dev/null || echo "N/A")
        dedup_count=$(samtools view -c "${sample_dir}/MQ25/dedup/s_all_${sample_name}_S_ChrM_MQ25_dedup.bam" 2>/dev/null || echo "N/A")

        echo "${sample_name},${chrm_count},${mq25_count},${dedup_count}" >> "$output_file"
        echo "Processed: $sample_name"
    done

    echo "Results for ChrM and subdirectories have been written to $output_file"
}

# Main execution
process_bam_files
process_chrm_directories

echo "Processing complete. Check the CSV files in ${base_dir}"