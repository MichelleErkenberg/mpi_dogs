import sys
import csv
from Bio import AlignIO

def process_alignment(aln_file, output_file, ref_id, file_format):
    # Read the alignment
    try:
        alignment = AlignIO.read(aln_file, file_format)
    except Exception as e:
        print(f"Error reading alignment file: {e}")
        sys.exit(1)

    # Create a dictionary for the sequences
    sequences = {record.id: str(record.seq) for record in alignment}

    # Check if the reference ID is present in the alignment
    if ref_id not in sequences:
        print(f"Error: Reference sequence with ID '{ref_id}' not found in the alignment.")
        sys.exit(1)

    # Get the reference sequence from the alignment
    ref_seq = sequences[ref_id]

    # Open the output file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write the header with the reference ID in parentheses
        header = ['Position', f'Reference ({ref_id})'] + [seq_id for seq_id in sequences if seq_id != ref_id]
        writer.writerow(header)

        # Write the data
        for i, ref_base in enumerate(ref_seq, start=1):
            row = [i, ref_base]
            for seq_id in sequences:
                if seq_id != ref_id:
                    row.append(sequences[seq_id][i-1])
            writer.writerow(row)

    print(f"Output written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <alignment_file> <output_file> <reference_sequence_id> <format>")
        print("Supported formats: fasta, clustal, phylip")
        sys.exit(1)

    aln_file = sys.argv[1]
    output_file = sys.argv[2]
    ref_id = sys.argv[3]
    file_format = sys.argv[4]

    process_alignment(aln_file, output_file, ref_id, file_format)
