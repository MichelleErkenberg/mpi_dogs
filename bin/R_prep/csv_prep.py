import pandas as pd
import sys
import re
import glob

def sort_sample_names(name):
    match = re.search(r'sample_(\d+)', name)
    return int(match.group(1)) if match else 0

def process_csv(input_pattern, output_file):
    # Get all CSV files matching the pattern
    input_files = glob.glob(input_pattern)
    
    all_results = []
    
    for file in input_files:
        # Extract the name from the file (assuming format all_env_Name.csv)
        env_name = re.search(r'all_env_(.+)\.csv', file).group(1)
        
        # Read the CSV file
        df = pd.read_csv(file)
        
        # Group the data by Sample and sum Matches and Total Reads
        grouped = df.groupby('Sample').agg({'Matches': 'sum', 'Total Reads': 'sum'})
        
        # Calculate the ratio of sums
        results = grouped['Matches'] / grouped['Total Reads']
        
        # Create a DataFrame with the results
        results_df = pd.DataFrame({env_name: results})
        
        all_results.append(results_df)
    
    # Concatenate all results
    final_df = pd.concat(all_results, axis=1)
    
    # Sort the index (Sample names)
    final_df.sort_index(key=lambda x: x.map(sort_sample_names), inplace=True)
    
    # Rename the index to 'Sample'
    final_df.index.name = 'Sample'
    
    # Save the results to the specified output CSV file
    final_df.to_csv(output_file)
    
    print(f"Results have been saved to {output_file}")

# Check if correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python script.py <input_pattern> <output_file>")
    sys.exit(1)

# Get input pattern and output file name from command line arguments
input_pattern = sys.argv[1]
output_file = sys.argv[2]

# Process the CSV files
process_csv(input_pattern, output_file)
