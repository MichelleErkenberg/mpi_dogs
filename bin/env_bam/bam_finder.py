import pysam
import csv
import argparse
import os

# Step 1: Set up argument parsing
parser = argparse.ArgumentParser(description="Search BAM files using positions from a CSV file.")
parser.add_argument('input_csv', type=str, help='The input CSV file containing positions and expected nucleotides.')
parser.add_argument('bam_file', type=str, help='The BAM file to search through.')
parser.add_argument('output_prefix', type=str, help='The prefix for the output CSV file.')
parser.add_argument('query_column', type=str, help='The name of the column to query for expected nucleotides.')

args = parser.parse_args()

# Step 2: Read positions and expected nucleotides from the input CSV file
positions = []
with open(args.input_csv, mode='r') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        # Append a tuple of (position, expected nucleotide) to the list
        position = int(row['Position'])
        expected_nucleotide = row[args.query_column]  # Use specified column for expected nucleotide
        positions.append((position, expected_nucleotide))

# Step 3: Open the BAM file for reading
bam_file = pysam.AlignmentFile(args.bam_file, "rb")
results = []

# Step 4: Iterate through each position and expected nucleotide
for position, expected_nucleotide in positions:
    # Perform pileup at the specified position (adjust chromosome name as needed)
    pileup_column = bam_file.pileup('chr1', position, position + 1)
    
    match_count = 0  # Count of matches with expected nucleotide
    total_count = 0  # Total reads at this position
    nucleotides_at_position = set()  # Set to store unique nucleotides found

    # Step 5: Analyze each pileup column
    for pileup in pileup_column:
        total_count += pileup.n  # Increment total read count
        
        for pileup_read in pileup.pileups:
            # Check if the base is not deleted or skipped
            if not pileup_read.is_del and not pileup_read.is_refskip:
                # Get the nucleotide from the query sequence
                nucleotide = pileup_read.alignment.query_sequence[pileup_read.query_position]
                nucleotides_at_position.add(nucleotide)  # Add nucleotide to the set
                
                # Check if it matches the expected nucleotide
                if nucleotide == expected_nucleotide:
                    match_count += 1

    # Append results for this position to the results list
    results.append({
        'Position': position,
        'Matches': match_count,
        'Total Reads': total_count,
        'Expected Nucleotide': expected_nucleotide,
        'Nucleotides Found': ', '.join(nucleotides_at_position)
    })

# Close the BAM file after processing
bam_file.close()

# Step 6: Write results to a new CSV file specified by the user
output_filename = f"{args.output_prefix}_{args.query_column}.csv"
with open(output_filename, mode='w', newline='') as outfile:
    fieldnames = ['Position', 'Matches', 'Total Reads', 'Expected Nucleotide', 'Nucleotides Found']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)

    writer.writeheader()  # Write header row
    for result in results:
        writer.writerow(result)  # Write each result as a row in the CSV file

print("Processing complete. Results saved to", output_filename)