import pysam
import csv
import argparse
import glob
import os

# Step 1: Set up argument parsing
parser = argparse.ArgumentParser(description="Search multiple BAM files using positions from a CSV file.")
parser.add_argument('input_csv', type=str, help='The input CSV file containing positions and expected nucleotides.')
parser.add_argument('bam_directory', type=str, help='The directory containing BAM files.')
parser.add_argument('output_file', type=str, help='The output CSV file.')
parser.add_argument('query_column', type=str, help='The name of the column to query for expected nucleotides.')

args = parser.parse_args()

# Step 2: Read positions and expected nucleotides from the input CSV file
positions = []
with open(args.input_csv, mode='r') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        position = int(row['Position'])
        expected_nucleotide = row[args.query_column].upper()
        positions.append((position, expected_nucleotide))

# Step 3: Get all BAM files in the specified directory
bam_files = glob.glob(os.path.join(args.bam_directory, "*.bam"))

# Step 4: Specify mitochondrial chromosome name
chromosome = "NC_002008.4"

# Step 5: Process each BAM file and write results to the output CSV
with open(args.output_file, mode='w', newline='') as outfile:
    fieldnames = ['Sample', 'Position', 'Matches', 'Total Reads', 'Expected Nucleotide', 'Nucleotides Found']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()

    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).split('.')[0]  # Extract sample name from BAM filename
        print(f"Processing {sample_name}...")

        samfile = pysam.AlignmentFile(bam_file, "rb")
        results = []

        for position, expected_nucleotide in positions:
            pileup_column = samfile.pileup(chromosome, position, position + 1)
            
            match_count = 0
            total_count = 0
            nucleotides_at_position = set()

            for pileup in pileup_column:
                total_count += pileup.n
                
                for pileup_read in pileup.pileups:
                    if not pileup_read.is_del and not pileup_read.is_refskip:
                        nucleotide = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                        nucleotides_at_position.add(nucleotide)
                        
                        if nucleotide == expected_nucleotide:
                            match_count += 1

            results.append({
                'Sample': sample_name,
                'Position': position,
                'Matches': match_count,
                'Total Reads': total_count,
                'Expected Nucleotide': expected_nucleotide,
                'Nucleotides Found': ', '.join(nucleotides_at_position)
            })

        samfile.close()

        # Write results for this sample
        for result in results:
            writer.writerow(result)
        
        # Add an empty row between samples
        writer.writerow({})

print(f"Processing complete. Results saved to {args.output_file}")
