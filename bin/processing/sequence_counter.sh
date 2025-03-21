#!/bin/bash

# Define the base directory
base_dir="${BASE_PATH}/data/dog_samples/processing"

# Function to extract the sample name from the filename
extract_sample_name() {
    local filename="$1"
    echo "$filename" | sed -n 's/.*s_all_\(.*\)_S.*/\1/p'
}

# Function to index BAM files in ChrM and subdirectories if index doesn't exist
index_bam() {
    local bam_file="$1"
    if [ ! -f "${bam_file}.bai" ]; then
        echo "Indexing ${bam_file}..."
        samtools index "$bam_file"
    fi
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
        chrm_file="${sample_dir}/s_all_${sample_name}_S_ChrM.bam"
        mq25_file="${sample_dir}/MQ25/s_all_${sample_name}_S_ChrM_MQ25.bam"
        dedup_file="${sample_dir}/MQ25/dedup/s_all_${sample_name}_S_ChrM_MQ25_dedup.bam"

        # Index BAM files only if necessary
        [ -f "$chrm_file" ] && index_bam "$chrm_file"
        [ -f "$mq25_file" ] && index_bam "$mq25_file"
        [ -f "$dedup_file" ] && index_bam "$dedup_file"

        chrm_count=$([ -f "$chrm_file" ] && samtools view -c "$chrm_file" || echo "N/A")
        mq25_count=$([ -f "$mq25_file" ] && samtools view -c "$mq25_file" || echo "N/A")
        dedup_count=$([ -f "$dedup_file" ] && samtools view -c "$dedup_file" || echo "N/A")

        echo "${sample_name},${chrm_count},${mq25_count},${dedup_count}" >> "$output_file"
        echo "Processed: $sample_name"
    done

    echo "Results for ChrM and subdirectories have been written to $output_file"
}

# Main execution
bam_file="${base_dir}/bam_files_sequence_counts.csv"
chrm_file="${base_dir}/chrm_sequence_counts.csv"

# Check and process bam_files
if [[ -f "$bam_file" ]]; then
    read -p "bam_files_sequence_counts.csv already exists. Repeat processing? (y/n): " choice
    if [[ $choice == "y" ]]; then
        process_bam_files
    else
        echo "Skipping bam_files processing."
    fi
else
    echo "Processing bam_files..."
    process_bam_files
fi

# Check and process chrm_files
if [[ -f "$chrm_file" ]]; then
    read -p "chrm_sequence_counts.csv already exists. Repeat processing? (y/n): " choice
    if [[ $choice == "y" ]]; then
        process_chrm_directories
    else
        echo "Skipping chrm_files processing."
    fi
else
    echo "Processing chrm_files..."
    process_chrm_directories
fi

echo "Processing complete. Check the CSV files in ${base_dir}"
