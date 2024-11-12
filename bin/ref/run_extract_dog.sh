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

# 2nd extraction for office 2 with 3 dogs (Cami as an outgroup) 
bash "$extract_script" "$FILE" "$OUTDIR/office_2/3dogs.csv" "Lily" "ThorA" "Cami"

#Anda as a substitute for ThorA
mkdir -p "$OUTDIR/office_2/Anda"
# 3rd extraction for office 2 with 2 dogs, using Anda as a substitute for ThorA
bash "$extract_script" "$FILE" "$OUTDIR/office_2/Anda/2dogs_Anda.as.ThorA.csv" "Lily" "Anda"

# 4th extraction for office 2 with 3 dogs (Cami as an outgroup) using Anda as a substitute for ThorA
bash "$extract_script" "$FILE" "$OUTDIR/office_2/Anda/3dogs_Anda.as.ThorA.csv" "Lily" "Anda" "Cami"
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
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_1/5dogs.csv" "$OUTDIR/office_1/5dogs.Cami.csv" "Cami"

# comparing all the dogs in office 2 for 2 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/2dogs.csv" "$OUTDIR/office_2/2dogs.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/2dogs.csv" "$OUTDIR/office_2/2dogs.ThorA.csv" "ThorA"

# comparing all the dogs in office 2 for 3 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/3dogs.csv" "$OUTDIR/office_2/3dogs.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/3dogs.csv" "$OUTDIR/office_2/3dogs.ThorA.csv" "ThorA"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/3dogs.csv" "$OUTDIR/office_2/3dogs.Cami.csv" "Cami"

#Anda as a substitute for ThorA
# comparing all the dogs in office 2 for 2 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/Anda/2dogs_Anda.as.ThorA.csv" "$OUTDIR/office_2/Anda/2dogs_Anda.as.ThorA.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/Anda/2dogs_Anda.as.ThorA.csv" "$OUTDIR/office_2/Anda/2dogs_Anda.as.ThorA.Anda.csv" "Anda"

# comparing all the dogs in office 2 for 3 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/Anda/3dogs_Anda.as.ThorA.csv" "$OUTDIR/office_2/Anda/3dogs_Anda.as.ThorA.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/Anda/3dogs_Anda.as.ThorA.csv" "$OUTDIR/office_2/Anda/3dogs_Anda.as.ThorA.Anda.csv" "Anda"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_2/Anda/3dogs_Anda.as.ThorA.csv" "$OUTDIR/office_2/Anda/3dogs_Anda.as.ThorA.Cami.csv" "Cami"
