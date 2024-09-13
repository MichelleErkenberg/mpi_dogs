#!/bin/bash

# Name of MSA file
msa_file="combined_cutoff.fasta"

# Name of the output CSV file
output_csv="dog_sequence_differences.csv"

# Function to extract the desired part of the dog name
extract_dog_name() {
    echo "$1" | sed -n 's/^>*s_all_\(.*\)_S8.*/\1/p'
}

# Read sequences and names from the MSA file
readarray -t full_sequences < "$msa_file"

# Extract names and sequences
names=()
sequences=()
current_sequence=""

for line in "${full_sequences[@]}"; do
    if [[ $line == ">"* ]]; then
        if [ ! -z "$current_sequence" ]; then
            sequences+=("$current_sequence")
            current_sequence=""
        fi
        full_name="${line}"
        short_name=$(extract_dog_name "$full_name")
        names+=("$short_name")
    else
        current_sequence+="$line"
    fi
done

# Add the last sequence
if [ ! -z "$current_sequence" ]; then
    sequences+=("$current_sequence")
fi

# Ensure we have the same number of names and sequences
if [ ${#names[@]} -ne ${#sequences[@]} ]; then
    echo "Error: Mismatch between number of names and sequences"
    exit 1
fi

# Initialize variables
seq_length=${#sequences[0]}
diff_count=0
diff_positions=()
diff_details=()

# Compare each position
for ((i=0; i<seq_length; i++)); do
    chars=()
    for seq in "${sequences[@]}"; do
        chars+=("${seq:$i:1}")
    done

    # Check for differences, ignore gaps (-)
    unique_chars=$(printf "%s\n" "${chars[@]}" | grep -v "-" | sort -u)
    if [ $(echo "$unique_chars" | wc -l) -gt 1 ]; then
        ((diff_count++))
        diff_positions+=($((i+1)))
        
        # Collect details about the differences
        detail="$((i+1))"
        for ((j=0; j<${#names[@]}; j++)); do
            detail+=",${chars[j]}"
        done
        diff_details+=("$detail")
    fi
done

# Output results to CSV file
echo "Position,${names[*]}" > "$output_csv"
for detail in "${diff_details[@]}"; do
    echo "$detail" >> "$output_csv"
done

# Output summary to console
echo "Number of different positions: $diff_count"
echo "Positions with differences: ${diff_positions[*]}"
echo "Detailed results have been saved to $output_csv"

# Visualization of differences
echo -e "\nVisualization of differences:"
for ((i=0; i<${#names[@]}; i++)); do
    echo -n "${names[i]}: "
    for ((j=0; j<seq_length; j++)); do
        char="${sequences[i]:j:1}"
        if [[ " ${diff_positions[*]} " =~ " $((j+1)) " ]]; then
            echo -ne "\033[31m$char\033[0m"
        else
            echo -n "$char"
        fi
    done
    echo
done

