# Overview
This repository contains data and scripts of the manuscript "Genetic variation reveals Individual-level climate tracking across the full annual cycle of a migratory bird"

# Description of datasets
Breeding_96snps(.loc,.ped,.map) - *PED formated genotypes from 96 assayed SNPs for all breeding individuals (including RAD sequenced*
Unknowns_96snps(.loc,.ped,.map) - *PED formated genotypes from 96 assayed SNPs for wintering and migrating birds*
Breeding_157snps.str - *STRUCTURE input with all breeding individuals assayed at 157 snps



# Description of scripts
Format.R - *filters and formats genotypes from Fluidigm assays for downstream analysis*. 
MakeMapStacks.R - *takes Q-matrix from STRUCTURE and creates spatial layers*. 
NewMap.R - *Plots layers from MakeMapStacks.R*. 
SNPrelate.R - *performs PCA*. 
Rubias.R - *performs population assignments in rubias*. 
OriGen.R - *creates probability surfaces for wintering individuals*. 
Morphology.R - *test for morphological correlates of precipitation with the Wiedenfeld 1991 data*. 
ClimateMatchingAnalysis.R - *test for correlation between breeding and wintering climate*. 

