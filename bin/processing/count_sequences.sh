##!/bin/bash


cd "$BASE_PATH/data/dog_samples"

# Define the base output directory
base_dir="./dog_samples"

# Function to process BAM files in a directory
process_directory() {
    local dir="$1"
    local output_file="$dir/sequence_counts.csv"

    # Check if there are any BAM files in this directory
    if ! ls "$dir"/*.bam 1> /dev/null 2>&1; then
        echo "No BAM files found in $dir"
        return
    fi

    # Write the header to the CSV file
    echo "Filename,SequenceCount" > "$output_file"

    # Iterate over all BAM files in the current directory
    for bam_file in "$dir"/*.bam
    do
        # Count the number of sequences in the BAM file
        sequence_count=$(samtools view -c "$bam_file")
        
        # Get the filename without the path
        filename=$(basename "$bam_file")

        # Print the filename and the number of sequences to the console
        echo "File: $filename - Number of sequences: $sequence_count"
        
        # Append the data to the CSV file
        echo "$filename,$sequence_count" >> "$output_file"
    done

    echo "Results for $dir have been written to $output_file"
}

# Use find to get all directories, including the base directory and all subdirectories
find "$base_dir" -type d | while read -r dir
do
    process_directory "$dir"
done

