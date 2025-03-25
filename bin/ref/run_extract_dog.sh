#!/bin/bash

# Path to the extraction script and reference coordinate file
extract_script="$BASE_PATH/bin/ref/extract.sh"
FILE="$BASE_PATH/data/dog_samples/ref/ref_coordinates.csv"
OUTDIR="$BASE_PATH/data/dog_samples/ref/"

while true; do

read -p "Use all dogs for office 1 and 2 OR filtering data to exclude dogs (raw/exclude)?: " x

if [[ "$x" == "raw" ]]; then
# creating directories for dog office 1 and 2	
mkdir -p "$OUTDIR/office_container"  #office 1
mkdir -p "$OUTDIR/office_thorA.lily" #office 2

#the references position is always extracted 
# First extraction for office 1 and 4 dogs
bash "$extract_script" "$FILE" "$OUTDIR/office_container/4dogs.csv" "Heidi" "Vito" "Urza" "Fritzy"
#second extraction for office 1 and 5 dogs (Cami as an outgroup)
bash "$extract_script" "$FILE" "$OUTDIR/office_container/5dogs.csv" "Heidi" "Vito" "Urza" "Fritzy" "Cami"
echo "extractions for office 1 finished"

# first extraction for office 2 with 2 dogs
bash "$extract_script" "$FILE" "$OUTDIR/office_thorA.lily/2dogs.csv" "Lily" "ThorA"
# second extraction for office 2 with 3 dogs (Cami as an outgroup) 
bash "$extract_script" "$FILE" "$OUTDIR/office_thorA.lily/3dogs.csv" "Lily" "ThorA" "Cami"
echo "extractions for office 1 finished"


#---------------second step - comparing dogs --------------------------------

echo "Continue to compare each dog office dog against each other."
#--------------office 1-------------
# Comparing all the dogs in office 1 for 4 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_container/4dogs.csv" "$OUTDIR/office_container/4dogs.Heidi.csv" "Heidi"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_container/4dogs.csv" "$OUTDIR/office_container/4dogs.Vito.csv" "Vito"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_container/4dogs.csv" "$OUTDIR/office_container/4dogs.Urza.csv" "Urza"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_container/4dogs.csv" "$OUTDIR/office_container/4dogs.Fritzy.csv" "Fritzy"

# comparing all the dogs in office 1 for 5 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_container/5dogs.csv" "$OUTDIR/office_container/5dogs.Heidi.csv" "Heidi"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_container/5dogs.csv" "$OUTDIR/office_container/5dogs.Vito.csv" "Vito"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_container/5dogs.csv" "$OUTDIR/office_container/5dogs.Urza.csv" "Urza"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_container/5dogs.csv" "$OUTDIR/office_container/5dogs.Fritzy.csv" "Fritzy"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_container/5dogs.csv" "$OUTDIR/office_container/5dogs.Cami.csv" "Cami"

# -----------------------office 2----------------------
# comparing all the dogs in office 2 for 2 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_thorA.lily/2dogs.csv" "$OUTDIR/office_thorA.lily/2dogs.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_thorA.lily/2dogs.csv" "$OUTDIR/office_thorA.lily/2dogs.ThorA.csv" "ThorA"

# comparing all the dogs in office 2 for 3 dogs
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_thorA.lily/3dogs.csv" "$OUTDIR/office_thorA.lily/3dogs.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_thorA.lily/3dogs.csv" "$OUTDIR/office_thorA.lily/3dogs.ThorA.csv" "ThorA"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/office_thorA.lily/3dogs.csv" "$OUTDIR/office_thorA.lily/3dogs.Cami.csv" "Cami"


#--------------- excluding dogs from the data ----------------

elif [[ "$x" == "exclude" ]]; then
	read -p "please decide whether to keep ThorA or Anda and ThorB or Cami: " a b 

#extracting all dogs (always decide between closely related once)
mkdir -p "$OUTDIR/all_dogs_with_${a}_${b}"
bash "$extract_script" "$FILE" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_with_${a}_${b}.csv" "$a" "$b" "Fritzy" "Heidi" "Urza" "Vito" "Lily" "Charlie" 

	read -p "also merge Lily into Anda/ThorA data (n/y)?: " l
	if [[ "l" == "y" ]]; then
mkdir -p "$OUTDIR/all_dogs_${a}_${b}_without_Lily"
bash "$extract_script" "$FILE" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}_${b}_without_Lily.csv" "$a" "$b" "Fritzy" "Heidi" "Urza" "Vito" "Charlie" 
echo "Extraction without Lily and $a and $b completed."
	else
echo "All extractions completed."
	fi
#-----------------comparing all dogs (minus closely related once)-------------------

python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_with_${a}_${b}.csv" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_${a}${b}.${a}.csv" "$a"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_with_${a}_${b}.csv" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_${a}${b}.${b}.csv" "$b"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_with_${a}_${b}.csv" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_${a}${b}.Heidi.csv" "Heidi"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_with_${a}_${b}.csv" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_${a}${b}.Fritzy.csv" "Fritzy"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_with_${a}_${b}.csv" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_${a}${b}.Vito.csv" "Vito"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_with_${a}_${b}.csv" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_${a}${b}.Urza.csv" "Urza"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_with_${a}_${b}.csv" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_${a}${b}.Lily.csv" "Lily"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_with_${a}_${b}.csv" "$OUTDIR/all_dogs_with_${a}_${b}/all_dogs_${a}${b}.Charlie.csv" "Charlie"

python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}_${b}_without_Lily.csv" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}${b}woL.${a}.csv" "$a"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}_${b}_without_Lily.csv" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}${b}woL.${b}.csv" "$b"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}_${b}_without_Lily.csv" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}${b}woL.Heidi.csv" "Heidi"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}_${b}_without_Lily.csv" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}${b}woL.Fritzy.csv" "Fritzy"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}_${b}_without_Lily.csv" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}${b}woL.Vito.csv" "Vito"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}_${b}_without_Lily.csv" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}${b}woL.Urza.csv" "Urza"
python3 "$BASE_PATH/bin/ref/diff_finder.py" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}_${b}_without_Lily.csv" "$OUTDIR/all_dogs_${a}_${b}_without_Lily/all_dogs_${a}${b}woL.Charlie.csv" "Charlie"
fi

read -p "Continue filtering (y/n)?: " q
	if [[ "$q" == "y" ]]; then
		continue
	else
		break
	fi	
done		 