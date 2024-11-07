#!/bin/bash

# Path to the extraction script and reference coordinate file
extract_script="$BASE_PATH/bin/ref/extract.sh"
FILE="$BASE_PATH/data/dog_samples/ref/ref_coordinates.csv"
OUTDIR="$BASE_PATH/data/dog_samples/ref/"

mkdir -p "$OUTDIR/office_1"
mkdir -p "$OUTDIR/office_2"

#the references position is always extracted 
# First extraction for office 1 and 4 dogs
bash "$extract_script" "$FILE" "$OUTDIR/office_1/4dogs.csv" "Heidi" "Vito" "Urza" "Fritzy"

#second extraction for office 1 and 5 dogs (Cami as an outgroup)
bash "$extract_script" "$FILE" "$OUTDIR/office_1/5dogs.csv" "Heidi" "Vito" "Urza" "Fritzy" "Cami"

# first extraction for office 2 with 2 dogs
bash "$extract_script" "$FILE" "$OUTDIR/office_2/2dogs.csv" "Lily" "ThorA"

# second extraction for office 2 with 3 dogs (Cami as an outgroup)
bash "$extract_script" "$FILE" "$OUTDIR/office_2/3dogs.csv" "Lily" "ThorA" "Cami"
echo "All extractions completed."

# Comparing all the dogs in office 1 for 4 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/4dogs.csv" "$OUTDIR/office_1/4dogs.Heidi.csv" "Heidi"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/4dogs.csv" "$OUTDIR/office_1/4dogs.Vito.csv" "Vito"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/4dogs.csv" "$OUTDIR/office_1/4dogs.Urza.csv" "Urza"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/4dogs.csv" "$OUTDIR/office_1/4dogs.Fritzy.csv" "Fritzy"

# comparing all the dogs in office 1 for 5 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/5dogs.csv" "$OUTDIR/office_1/5dogs.Heidi.csv" "Heidi"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/5dogs.csv" "$OUTDIR/office_1/5dogs.Vito.csv" "Vito"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/5dogs.csv" "$OUTDIR/office_1/5dogs.Urza.csv" "Urza"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/5dogs.csv" "$OUTDIR/office_1/5dogs.Fritzy.csv" "Fritzy"

# comparing all the dogs in office 2 for 2 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/2dogs.csv" "$OUTDIR/office_2/2dogs.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/2dogs.csv" "$OUTDIR/office_2/2dogs.ThorA.csv" "ThorA"

# comparing all the dogs in office 2 for 3 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/3dogs.csv" "$OUTDIR/office_2/3dogs.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/3dogs.csv" "$OUTDIR/office_2/3dogs.ThorA.csv" "ThorA"
