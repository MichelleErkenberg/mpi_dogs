import pysam
import csv
import argparse
import glob
import os
import re

# Function to sort sample names numerically
def sort_sample_names(name):
    match = re.search(r'sample_(\d+)', name)
    if match:
        return int(match.group(1))
    return 0

# Step 1: Set up argument parsing
parser = argparse.ArgumentParser(description="Search multiple BAM files using positions from a CSV file.")
parser.add_argument('input_csv', type=str, help='The input CSV file containing positions and expected nucleotides.')
parser.add_argument('bam_directory', type=str, help='The directory containing BAM files.')
parser.add_argument('output_file', type=str, help='The output CSV file.')
parser.add_argument('query_column', type=str, help='The name of the column to query for expected nucleotides.')
parser.add_argument('--chromosome', type=str, default=None, help='The chromosome name to use. If not provided, will use the first chromosome in each BAM file.')
parser.add_argument('--quality_threshold', type=int, default=20, help='The quality threshold for base calls. Default is 20.')

args = parser.parse_args()

# Step 2: Read positions and expected nucleotides from the input CSV file
positions = []
with open(args.input_csv, mode='r') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        position = int(row['Position'])
        expected_nucleotide = row[args.query_column].upper()
        positions.append((position, expected_nucleotide))

# Step 3: Get all BAM files in the specified directory and sort them
bam_files = glob.glob(os.path.join(args.bam_directory, "*.bam"))
bam_files.sort(key=lambda x: sort_sample_names(os.path.basename(x)))

# Step 4: Process each BAM file and write results to the output CSV
with open(args.output_file, mode='w', newline='') as outfile:
    fieldnames = ['Sample', 'Chromosome', 'Position', 'Matches', 'Total Reads', 'Expected Nucleotide', 'Nucleotides Found', 'Deletions', 'Insertions', 'N_count', 'Low_quality', 'Softclipped']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()

    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).split('.')[0]  # Extract sample name from BAM filename
        print(f"Processing {sample_name}...")

        samfile = pysam.AlignmentFile(bam_file, "rb")
        
        # Determine which chromosome to use
        if args.chromosome:
            chromosome = args.chromosome
        else:
            chromosome = next(iter(samfile.references))  # Use the first chromosome in the BAM file
        
        results = []

        for position, expected_nucleotide in positions:
            try:
                pileup_column = samfile.pileup(chromosome, position - 1, position)  # Adjust for 0-based indexing
                
                match_count = 0
                total_count = 0
                nucleotides_at_position = set()
                deletions = 0
                insertions = 0
                n_count = 0
                low_quality = 0
                softclipped = 0

                for pileup in pileup_column:
                    if pileup.pos == position - 1:  # Check if we're at the correct position
                        total_count = pileup.n
                        
                        for pileup_read in pileup.pileups:
                            if pileup_read.is_del:
                                deletions += 1
                            elif pileup_read.is_refskip:
                                continue
                            else:
                                query_pos = pileup_read.query_position
                                if query_pos is not None:
                                    nucleotide = pileup_read.alignment.query_sequence[query_pos].upper()
                                    nucleotides_at_position.add(nucleotide)
                                    
                                    if nucleotide == 'N':
                                        n_count += 1
                                    elif pileup_read.alignment.query_qualities[query_pos] < args.quality_threshold:
                                        low_quality += 1
                                    elif nucleotide == expected_nucleotide:
                                        match_count += 1
                                
                                if pileup_read.indel > 0:
                                    insertions += 1
                            
                            if pileup_read.alignment.cigartuples[0][0] == 4:  # Softclipped at start
                                softclipped += 1
                            if pileup_read.alignment.cigartuples[-1][0] == 4:  # Softclipped at end
                                softclipped += 1

                results.append({
                    'Sample': sample_name,
                    'Chromosome': chromosome,
                    'Position': position,
                    'Matches': match_count,
                    'Total Reads': total_count,
                    'Expected Nucleotide': expected_nucleotide,
                    'Nucleotides Found': ', '.join(nucleotides_at_position),
                    'Deletions': deletions,
                    'Insertions': insertions,
                    'N_count': n_count,
                    'Low_quality': low_quality,
                    'Softclipped': softclipped
                })
            except ValueError:
                # This will catch cases where the chromosome or position is not found in the BAM file
                results.append({
                    'Sample': sample_name,
                    'Chromosome': chromosome,
                    'Position': position,
                    'Matches': 0,
                    'Total Reads': 0,
                    'Expected Nucleotide': expected_nucleotide,
                    'Nucleotides Found': 'Position not found',
                    'Deletions': 0,
                    'Insertions': 0,
                    'N_count': 0,
                    'Low_quality': 0,
                    'Softclipped': 0
                })

        samfile.close()

        # Write results for this sample
        for result in results:
            writer.writerow(result)
        
        # Add an empty row between samples
        writer.writerow({})

print(f"Processing complete. Results saved to {args.output_file}")
