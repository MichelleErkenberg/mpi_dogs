#!/bin/bash

# Define the base directory
base_dir="${BASE_PATH}/data/dog_samples/processing"

# Function to extract the sample name from the filename
extract_sample_name() {
    local filename="$1"
    echo "$filename" | sed -n 's/.*s_all_\(.*\)_S.*/\1/p'
}

# Function to index BAM files in ChrM directory (only if .bai doesn't exist)
index_bam_files_in_chrm() {
    find "${base_dir}/ChrM" -type f -name "*ChrM*.bam" | while read -r bam_file; do
        if [ ! -f "${bam_file}.bai" ]; then
            samtools index "$bam_file"
        fi
    done
}

# Function to process ChrM BAM files
process_chrm_files() {
    local dir="$base_dir/ChrM"
    local output_file="$base_dir/chrm_sequence_counts.csv"

    # Check if there are any ChrM BAM files in this directory
    if ! ls "$dir"/*ChrM*.bam 1> /dev/null 2>&1; then
        echo "No ChrM BAM files found in $dir"
        return
    fi

    # Write the header to the CSV file
    echo "Sample,SequenceCount" > "$output_file"

    # Iterate over all ChrM BAM files in the directory
    for bam_file in "$dir"/*ChrM*.bam
    do
        # Count the number of sequences in the BAM file
        sequence_count=$(samtools view -c "$bam_file")
        
        # Extract the sample name
        sample_name=$(extract_sample_name "$(basename "$bam_file")")

        # Append the data to the CSV file
        echo "$sample_name,$sequence_count" >> "$output_file"
        
        echo "Processed: $sample_name - Count: $sequence_count"
    done

    echo "Results have been written to $output_file"
}

# Main execution
chrm_file="$base_dir/chrm_sequence_counts.csv"

# Check if the output file already exists
if [[ -f "$chrm_file" ]]; then
    read -p "chrm_sequence_counts.csv already exists. Repeat processing? (y/n): " choice
    if [[ $choice != "y" ]]; then
        echo "Skipping processing."
        exit 0
    fi
fi

# Index BAM files in ChrM directory before processing (only if necessary)
index_bam_files_in_chrm

# Process ChrM files
process_chrm_files

echo "Processing complete. Check the CSV file in ${base_dir}"