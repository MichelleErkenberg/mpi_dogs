from Bio import AlignIO
import csv
import re
import sys

# Check if both input and output filenames were provided as command-line arguments
if len(sys.argv) < 3:
    print("Usage: python dog_vari.py <input_alignment_file> <output_csv_file>")
    sys.exit(1)

# Get the input and output filenames from the command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the alignment
try:
    alignment = AlignIO.read(input_file, "fasta")
except FileNotFoundError:
    print(f"Error: File '{input_file}' not found.")
    sys.exit(1)
except Exception as e:
    print(f"Error reading file: {e}")
    sys.exit(1)

# Open CSV file for writing
with open(output_file, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    
    # Create header row with modified sequence names
    header = ["Position"]
    for record in alignment:
        if record.id.startswith("s_all_") and "_S" in record.id:
            match = re.search(r's_all_(.+)_S', record.id)
            if match:
                header.append(match.group(1))
            else:
                header.append(record.id)
        else:
            header.append(record.id)
    
    csvwriter.writerow(header)

    # Compare sequences
    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        unique_bases = set(column)
        
        # Check if there are exactly two unique bases in the column
        if len(unique_bases) == 2:
            # Check if the difference is only one base
            base_counts = {base: column.count(base) for base in unique_bases}
            if min(base_counts.values()) == 1:
                row = [i+1]  # Position
                for base in column:
                    row.append(base)
                csvwriter.writerow(row)

print(f"Processing complete. Check '{output_file}'.")
