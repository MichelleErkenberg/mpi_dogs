#!/bin/ bash

#script uses samtools to get for every bam file an bai file, important for further work

for INFILE in "$@"
do
	samtools index $INFILE
done
