#!/bin/bash

#creating and define new env_bam file as csv output file
mkdir -p "$BASE_PATH/data/dog_samples/env_bam"
env_bam="$BASE_PATH/data/dog_samples/env_bam"
export bam_file="$BASE_PATH/data/env_samples/quicksand.v2/out/Canidae/fixed/3-deduped/"
office_file="$BASE_PATH/data/dog_samples/ref/"

#creating *.bam.bai data for every bam to process
#bash env_bam/bam_to_bai.sh

#creating seperated file for our dogs to structure
#mkdir -p "$env_bam/Heidi"

#finding 
python3 env_bam/bam_finder_new.py "$office_file/office_1/5dogs.Heidi.csv" "$office_file/ref_coordinates.csv" "$bam_file" "$env_bam/all_env_Heidi.csv" "Heidi"
python3 env_bam/bam_finder_new.py "$office_file/office_1/5dogs.Fritzy.csv" "$office_file/ref_coordinates.csv" "$bam_file" "$env_bam/all_env_Fritzy.csv" "Fritzy"
python3 env_bam/bam_finder_new.py "$office_file/office_1/5dogs.Vito.csv" "$office_file/ref_coordinates.csv" "$bam_file" "$env_bam/all_env_Vito.csv" "Vito"
python3 env_bam/bam_finder_new.py "$office_file/office_1/5dogs.Urza.csv" "$office_file/ref_coordinates.csv" "$bam_file" "$env_bam/all_env_Urza.csv" "Urza"

python3 env_bam/bam_finder_new.py "$office_file/office_2/3dogs.Lily.csv" "$office_file/ref_coordinates.csv" "$bam_file" "$env_bam/all_env_Lily.csv" "Lily"
python3 env_bam/bam_finder_new.py "$office_file/office_2/3dogs.ThorA.csv" "$office_file/ref_coordinates.csv" "$bam_file" "$env_bam/all_env_ThorA.csv" "ThorA"
