library(data.table)
library(ggplot2)

# if (is_student)
#   setwd('//biofs09/scratch/praktikum2024_mol_anth/')
# else
#   setwd('/Volumes/scratch/praktikum2024_mol_anth/')

setwd("~/github/mpi_dogs/")


## Read in the sample sheet
dt.samples <- fread('data/dog_samples/R_prep/dog_env_samples_24_v1.txt', na.strings = c('-','NA',''))
setnames(dt.samples, 'dog', 'dogs_in_office')

## DIAGNOSTIC POSITIONS just the container dogs
#dt.dogs.container <- fread('data/container_dogs_env_pct_merged_data.tsv', na.strings = c('-','NA',''))

## DIAGNOSTIC POSITIONS all mpi dogs
#dt.dogs.all <- fread('classwork/schul01/clean01/target_dogs.aln.diag.frac.txt', na.strings = c('-','NA',''))
#setnames(dt.dogs.all, 'sample_1', 'sample_id')
#dt.dogs.all <- merge(dt.dogs.all, dt.samples, by='sample_id', all=T)

## DIAGNOSTIC POSITIONS all mpi dogs NO THOR
#dt.dogs.no_thor <- fread('classwork/schul01/clean01/target_dogs_no_Thors.aln.diag.frac.txt', na.strings = c('-','NA',''))
#setnames(dt.dogs.no_thor, 'sample_1', 'sample_id')
#dt.dogs.no_thor <- merge(dt.dogs.no_thor, dt.samples, by='sample_id', all=T)

## DIAGNOSTIC POSITIONS all mpi dogs NO ANDA or CAMI
#dt.dogs.no_a_c <- fread('classwork/schul01/clean01/target_dogs_no_Anda_Cami.aln.diag.frac.txt', na.strings = c('-','NA',''))
#setnames(dt.dogs.no_a_c, 'sample_1', 'sample_id')
#dt.dogs.no_a_c <- merge(dt.dogs.no_a_c, dt.samples, by='sample_id', all=T)

## DIAGNOSTIC POSITIONS all mpi dogs TARGETING ANDA VS CHARLIE
#dt.dogs.ac <- fread('classwork/schul01/clean01/target_dogs_only_Anda_Charlie.aln.diag.frac.txt', na.strings = c('-','NA',''))
#setnames(dt.dogs.ac, 'sample_1', 'sample_id')
#dt.dogs.ac <- merge(dt.dogs.ac, dt.samples, by='sample_id', all=T)

## TAXONOMIC GROUPS found in each sample
dt.tax <- fread('data/env_samples/quicksand.v2/final_report.tsv', na.strings = c('-','NA',''))
#setnames(dt.tax, 'dog', 'dogs_in_office')
##dog isn't part of the tsv file

## CLEAN TAXONOMIC GROUPS - make a table where every sample is listed four times, one for each Family
dt.cj_samples <- CJ(Family = c('Hominidae', 'Canidae', 'Felidae', 'Suidae'),
                    sample_id = unique(dt.tax$sample_id))
dt.cj_samples <- merge(dt.cj_samples, unique(dt.tax[, .(sample_id, ReadsRaw)]),
                       by = 'sample_id')

dt.tax.clean <- merge(dt.tax[, .(sample_id, ReadsDeduped, Family)], 
                      dt.cj_samples, 
                      allow.cartesian = T, 
                      all.y=T, by=c('Family', 'sample_id'))
dt.tax.clean[is.na(ReadsDeduped), ReadsDeduped := 0]
dt.tax.clean <- merge(dt.tax.clean, dt.samples, by='sample_id', all=T)
dt.tax.clean[is.na(y), y := 0]
dt.tax.clean[, ReadsDeduped.cap := ReadsDeduped]
dt.tax.clean[ReadsDeduped > 100000, ReadsDeduped.cap := 20000]




