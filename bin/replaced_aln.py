import argparse
import csv

def read_csv(filename):
    """
    Read the CSV file and organize data by sequence.
    
    The CSV file is expected to have positions in the first column
    and sequence data in subsequent columns.
    """
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        headers = next(reader)  # Read the header row
        # Create a dictionary for each sequence, excluding the 'Position' column
        data = {headers[i]: [] for i in range(1, len(headers))}
        for row in reader:
            position = int(row[0])
            # Iterate through each sequence's data, starting from column 1
            for i, value in enumerate(row[1:], start=1):
                if value != 'n':  # Ignore 'n' values as they represent no change
                    data[headers[i]].append((position, value))
    return data

def read_aln(filename):
    """Read the ALN file and return its contents as a list of lines."""
    with open(filename, 'r') as f:
        return f.readlines()

def write_aln(filename, lines):
    """Write the updated alignment to a new ALN file."""
    with open(filename, 'w') as f:
        f.writelines(lines)

def update_sequence(sequence, csv_data):
    """
    Update a single sequence with data from the CSV file.
    
    :param sequence: The original sequence string
    :param csv_data: List of (position, value) tuples for this sequence
    :return: Updated sequence string
    """
    sequence_list = list(sequence)
    for pos, value in csv_data:
        if 1 <= pos <= len(sequence_list):
            sequence_list[pos - 1] = value  # Adjust for 0-based indexing
    return ''.join(sequence_list)

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Update Alignment File with CSV Data")
    parser.add_argument("csv_file", help="Path to CSV file with replaced sequences")
    parser.add_argument("aln_file", help="Path to original ALN file")
    parser.add_argument("output_file", help="Path to new ALN file")
    args = parser.parse_args()

    # Read input files
    csv_data = read_csv(args.csv_file)
    aln_lines = read_aln(args.aln_file)

    # Process and update sequences in the ALN file
    updated_aln = []
    current_sequence = ""
    current_header = ""
    sequence_index = 0

    for line in aln_lines:
        if line.startswith('>'):
            # Process the previous sequence if it exists
            if current_sequence:
                if sequence_index < len(csv_data):
                    csv_key = list(csv_data.keys())[sequence_index]
                    updated_sequence = update_sequence(current_sequence, csv_data[csv_key])
                    # Split the updated sequence into 60-character lines
                    for i in range(0, len(updated_sequence), 60):
                        updated_aln.append(updated_sequence[i:i+60] + '\n')
                else:
                    # If no CSV data for this sequence, keep it unchanged
                    updated_aln.append(current_sequence)
                sequence_index += 1

            # Start processing a new sequence
            updated_aln.append(line)  # Add the header line
            current_header = line.strip()[1:]  # Store the header without '>'
            current_sequence = ""
        else:
            # Accumulate the sequence data
            current_sequence += line.strip()

    # Process the last sequence
    if current_sequence:
        if sequence_index < len(csv_data):
            csv_key = list(csv_data.keys())[sequence_index]
            updated_sequence = update_sequence(current_sequence, csv_data[csv_key])
            for i in range(0, len(updated_sequence), 60):
                updated_aln.append(updated_sequence[i:i+60] + '\n')
        else:
            updated_aln.append(current_sequence)

    # Save the updated alignment
    write_aln(args.output_file, updated_aln)

    print(f"Updated alignment has been saved to {args.output_file}.")

if __name__ == "__main__":
    main()
