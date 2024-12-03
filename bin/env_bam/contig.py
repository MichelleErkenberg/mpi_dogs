import pysam
import glob
import os

# Directory containing the BAM files
bam_directory = "../data/env_samples/Canidae"  

# Output file name
output_file = "contig_list.txt" 


# Find all BAM files in the specified directory
bam_files = glob.glob(os.path.join(bam_directory, "*.bam"))

# Open the output file in write mode
with open(output_file, 'w') as outfile:
    # Iterate through each BAM file and write the contig names to the file
    for bam_file in bam_files:
        # Write the name of the BAM file being processed
        outfile.write(f"Processing {os.path.basename(bam_file)}...\n")
        try:
            # Open the BAM file
            samfile = pysam.AlignmentFile(bam_file, "rb")
            # Get the list of contigs
            contigs = samfile.references
            outfile.write("Available contigs:\n")
            # Write each contig name to the file
            for ref in contigs:
                outfile.write(f"{ref}\n")
            # Close the BAM file
            samfile.close()
        except Exception as e:
            # If there's an error processing the BAM file, write the error message
            outfile.write(f"Error processing {bam_file}: {e}\n")
        # Add a blank line between BAM files for readability
        outfile.write("\n")

    # Write a completion message to the file
    outfile.write("Contig listing complete.\n")

# Print a confirmation message to the console
print(f"Contig listing complete. Results saved to {output_file}")
