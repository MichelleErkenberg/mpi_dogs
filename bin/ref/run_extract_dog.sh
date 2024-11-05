#!/bin/bash

# Path to the extraction script and reference coordinate file
extract_script="$BASE_PATH/bin/ref/extract.sh"
FILE="$BASE_PATH/data/dog_samples/ref/ref_coordinates.csv"
OUTDIR="$BASE_PATH/data/dog_samples/ref/"

mkdir -p "$OUTDIR/office_1"
mkdir -p "$OUTDIR/office_2"

#the references position get allways extracted 
# First extraction
bash "$extract_script" "$FILE" "$OUTDIR/office_1/4dogs.csv" "Heidi" "Vito" "Urza" "Fritzy"

# Second extraction
bash "$extract_script" "$FILE" "$OUTDIR/office_2/2dogs.csv" "Lily" "ThorA"

echo "All extractions completed."

# Comparing all the dogs in an office
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/4dogs.csv" "$OUTDIR/office_1/4dogs.Heidi.csv" "Heidi"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/4dogs.csv" "$OUTDIR/office_1/4dogs.Vito.csv" "Vito"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/4dogs.csv" "$OUTDIR/office_1/4dogs.Urza.csv" "Urza"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/4dogs.csv" "$OUTDIR/office_1/4dogs.Fritzy.csv" "Fritzy"

python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/2dogs.csv" "$OUTDIR/office_2/2dogs.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/2dogs.csv" "$OUTDIR/office_2/2dogs.ThorA.csv" "ThorA"
