#!/bin/bash

# Check if enough arguments are provided
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <input_file> <output_file> <column1> [column2] ..."
    exit 1
fi

# First two arguments are stored as input and output file names
input_file="$1"
output_file="$2"
# Shift the first two arguments out of the list 
shift 2

# Define the columns to extract; first column is fixed and get always extracted (reference coordinates); and additional arguments can be extracted (defined in the real run of the script)
columns_to_extract=("1" "$@")

# Function to extract columns from a CSV file, excluding lines with 'n'
extract_columns() {
    local columns=$(IFS=,; echo "${columns_to_extract[*]}")
    
    # Use grep to filter out lines containing 'n' (except for the header), then use csvcut to extract the specified columns
    {
        # Output the header
        head -n 1 "$input_file"
        # Output the filtered content
        tail -n +2 "$input_file" | grep -v 'n'
    } | csvcut -c "$columns" > "$output_file"
    
    echo "Extraction completed. Result saved in $output_file"
}

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file $input_file not found."
    exit 1
fi

# Perform the extraction
extract_columns
