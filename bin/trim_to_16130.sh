#!/bin/bash

# Set the working directory
work_dir="/mnt/expressions/michelle_erkenberg/github/mpi_dogs/data/dog_samples/Consensus/mask"
cd "$work_dir" || exit 1

# Set the desired length
desired_length=16130
# Set the number of bases per line (adjust this to match your original files)
bases_per_line=60

# create a 'trim' directory if it doesn't exist
mkdir -p trim

# Function to get the actual sequence name from the FASTA file
get_seq_name() {
    head -n 1 "$1" | sed 's/^>//' | awk '{print $1}' | sed 's/:[0-9-]*$//'
}

# Function to format FASTA sequence
format_fasta() {
    local input_file="$1"
    local output_file="$2"
    local seq_name="$3"
    
    # Write the sequence name
    echo ">$seq_name" > "$output_file"
    
    # Write the sequence data
    awk -v l=$bases_per_line 'NR>1 {
        for (i=1; i<=length($0); i+=l)
            print substr($0, i, l)
    }' "$input_file" >> "$output_file"
}

# For each .fas file in the current directory (mask)
for file in *.fas; do
    # Check if the file exists
    if [ -f "$file" ]; then
        # Extract the filename without extension
        base_name="${file%.fas}"

        # Get the actual sequence name from the FASTA file
        seq_name=$(get_seq_name "$file")
 
        # Create a temporary .bed file
        echo -e "${seq_name}\t0\t${desired_length}" > "${base_name}.bed"

        # Use bedtools to trim the sequence
        bedtools getfasta -fi "$file" -bed "${base_name}.bed" -fo "${work_dir}/trim/${base_name}.temp.fas"

        # Format the trimmed FASTA file
        format_fasta "${work_dir}/trim/${base_name}.temp.fas" "${work_dir}/trim/${base_name}.trimmed.fas" "$seq_name"

        # Remove the temporary files
        rm "${base_name}.bed"
        rm "${work_dir}/trim/${base_name}.temp.fas"

        echo "File $file has been trimmed to ${desired_length} bases and saved as ${work_dir}/trim/${base_name}.trimmed.fas"
    fi
done

echo "All .fas files have been processed. Trimmed files are in the '${work_dir}/trim' directory."
