#!/bin/bash

# Path to the extraction script and reference coordinate file
extract_script="./ref_test.sh"

#the reference and position in the reference coordinated are allways extracted 
# First extraction
"$extract_script" "data/ref_coordinates.csv" "data/3dogs.csv" "Anda" "Cami" 

echo "All extractions completed."

python3 diff_finder.py "data/3dogs.csv" "data/3dogs_diff_Anda.csv" "Anda"
