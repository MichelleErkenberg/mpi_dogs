#!/bin/bash

#change to directory with the *.bam files
cd ../data/dog_samples

for file in s_all_*_S*.bam; do
    # extracts the dog name
    newname=$(echo $file | sed 's/s_all_\(.*\)_S.*/\1/')
    
    # does that for all *.bam files
    newname="${newname}.bam"
    
    # renames the long file name in just our dogs file name
    mv "$file" "$newname"
    
    echo "renamed: $file -> $newname"
done
