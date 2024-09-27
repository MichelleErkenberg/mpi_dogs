##!/bin/bash

#change to correct dir

cd ../data/dog_samples

for file in s_all_*_S*.bam; do
    # extract the dog name
    newname=$(echo $file | sed 's/s_all_\(.*\)_S.*/\1/')
    
    # does this for all *.bam files
    newname="${newname}.bam"
    
    # renames the file
    mv "$file" "$newname"
    
    echo "Umbenannt: $file -> $newname"
done
