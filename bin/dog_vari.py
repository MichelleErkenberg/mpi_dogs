from Bio import AlignIO
import csv
import re

# Read the alignment
alignment = AlignIO.read("related_seq.trimmed.aln", "fasta")

# Open CSV file for writing
with open("sequence_differences.csv", "w", newline="") as csvfile:
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
        if len(set(column)) > 1:
            row = [i+1]  # Position
            for base in column:
                row.append(base)
            csvwriter.writerow(row)

print("Processing complete. Check 'sequence_differences.csv'.")
