#!/bin/bash

# Set the desired length
desired_length=16130

# create a 'trim' directory if it doesn't exist
mkdir -p trim

# For each .fas file in the current directory (mask)
for file in [A-Z]*.fas; do
    # Check if the file exists
    if [ -f "$file" ]; then
        # Extract the filename without extension
        base_name="${file%.fas}"

	  # Get the actual sequence name from the FASTA file
        seq_name=$(head -n 1 "$file" | sed 's/^>//')       
 
        # Create a temporary .bed file
        echo -e "${seq_name}\t0\t${desired_length}" > "${base_name}.bed"
        
        # Use bedtools to trim the sequence
        bedtools getfasta -fi "$file" -bed "${base_name}.bed" -fo "trim/${base_name}.trimmed.fas"
        
        # Remove the temporary .bed file
        rm "${base_name}.bed"
        
        echo "File $file has been trimmed to ${desired_length} bases and saved as trim/${base_name}.trimmed.fas"
    fi
done

echo "All .fas files have been processed. Trimmed files are in the 'trim' directory."
