# MPI Dog Office Project

This project was conducted with eight canines residing in their respective owners' offices at the MPI-EVA in Leipzig.The objective of this project was to identify genetic variations among the canines using DNA analysis. Leveraging this genomic information, the relationship among the canines was ascertained. With the assistance of reference canines, unsequenced information was annotated the subjects. Furthermore, samples were collected from various locations at the MPI-EVA, including office floors, hallways, elevators, and laboratories. The extracted DNA was then used to identify the locations where each dog was present. 

## One script to rule them all

For convenience, all the necessary steps for processing the data are consolidated into a single script:
```
bash /bin/run_all.sh
```

To use this in your environment, you need to change the BASE_PATH to:
```
export BASE_PATH="/Path/to/your/file" 
```



## Step-By-Step explanation of the script

### 1. Processing

Here it is nessesary to change the *.bam file to your *.bam file directory.
In a first step the data was processed. Therefore the raw data was filtered for the mitochondrial DNA (mtDNA). Furthermore a map quality of 25 % was required, and only deduplicated data was allowed.
```
bash processing/filter.sh
```

Also there are two scripts to count the sequences. One to count the amount of sequences in all chromosoms and one just for the amount of sequences that are mtDNA. This step is optional.


### 2. Consensus 

Using the processed data to call a consensus sequence for each dog is done by using Matthias' perl script and a canines references FASTA.
```
bash consensus/consensus.sh	
```

### 3. MSA

In a first step to create a multiple sequence alignment (MSA), a masking step was performed. This resulted in a renamed copy of the consensus sequence for all dogs in a separate file. In addition, the undetermined data were deleted and the files are renamed to a shorter form.
```
bash msa/mask.sh	
```

For further processing, all dog sequences go through a cutoff process. Therefore, all sequences were trimmed after a TTTTAGG/AAG part at the end of their DNA. 
```
bash msa/cutoff_seq.sh	
```

In a final step, the sequences of all our dogs (combined), previously published ones (combined_pub), and a reference dog (combined_pub.ref) were aligned using MAFFT. This resulted in three different .fasta and .aln files, shown in brackets.
```
bash msa/msa.sh 	
```

### 4. Phylogenetic tree

Using the alignment files, a phylogenetic tree was generated for the MPI dogs only, and for the MPI dogs, the previously published dogs, and the reference dog. Therefore FastTree was used.
```
bash tree/tree.sh
```

### 5. Detecting differences between the MPI dogs

The information from the phylogenetic tree was used to find closely related dogs for each MPI dog. This information was then manually added to the /bin/diff/related.txt file. Afterward the script can be executed.
```
bash diff/diff.sh
```

In a first step, the script trims the sequences of the file combined_pub.ref.aln to a length of 16138 bp. From these trimmed sequences, the sequences of the MPI dogs and their closely related dogs are extracted. The information is then transferred to a csv file. In this csv file, the n's of the MPI dogs are replaced by the bases of their closely related dogs, if they are the same. The replaced sequences of the MPI dogs were then written into a fasta file together with the sequence information of the reference dog. A new alignment was then performed.

### 6. Genome Coordinates

In order to compare the sequence of the MPI dogs with the sequence information from the environmental data, genomic coordinates are required. Therefore, the reference dog genome is used to assign coordinates to all MPI dog sequences.
```
python3 ref/ref_coor.py "$BASE_PATH/data/dog_samples/diff/replaced_seq.related_n.mpi_dogs.added_ref.aln" "$BASE_PATH/data/dog_samples/ref/ref_coordinates.csv" NC_002008.4 fasta
```

### 7. Detecting private positions

As we assume that MPI dogs are most likely to be found in their owner's office, we focused mainly on the differentiation between dogs in the same office. The following scripts only allow to compare dogs living in the same office to detect private positions. The offices of the dogs were added manually in the first script.
```
bash ref/run_extract_dog.sh 
bash env_bam/run_bam.sh
```

### 8. R preparation

In order to visualize the collected data in R, some processing was required. Therefore, a csv file was created containing the average radio for each MPI dog in each sample.
```
python3 R_prep/csv_prep.py "$BASE_PATH/data/dog_samples/env_bam/all_env_*.csv" "$BASE_PATH/data/dog_samples/R_prep/R_prep_sample_vs_dog.csv"
```

## Python packages

* pandas
* csv
* argparse
* sys
* re
* ete3
* pysam
* glob
* os

## R

R was used for data visualization. For further data analysis, we decided to group genetically similar dogs together. This includes Anda, Thor A and Lily as well as Thor B and Cami.

### Locations of the dogs

Samples were collected from different locations. To highlight the amount of dog DNA at each location, the following R script can be used
```
R/dog_catagorized_readsdeduped.R
```
As the data contains samples with no data points or low amounts of DNA, these samples need to be sorted out. For SNP validation, the position of the SNPs must be covered by at least 10 sequences, regardless of whether the SNP occurs in them or not. There must also be at least 50 ReadsDeduped in the raw data. 
After this filtering, different locations where the dogs might be present are examined. 

Samples that are filtered out are highlighted with an asterisk. The data shows that each dog is mainly present in its owner's office. Based on this knowledge, further research was conducted.

### Headmaps

With the knowledge that every dog's DNA is most likely to occur in there owners office, we started to investigate the two main dog offices. Those include for office 1: Heidi, Vito, Fritzy and Urza, and for office 2: Lily and Thor A. As Lily and Thor A were grouped due to there genetic similarities in the first place, we used a different filtered dataset. The amount of DNA belonging to each individual dog was visualized as a heatmap.

### Different species in dog and non-dog locations

As the ground samples were taken from various locations, they might as well contain DNA from other species. Therefore, the presence of human, Felican and Suidae DNA was also checked.

