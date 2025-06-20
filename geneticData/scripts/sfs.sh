#!/bin/bash

## This script is meant to calculate the site frequency spectrum of a merged vcf file

# Load vcftools
module load vcftools-0.1.14-gcc-11.2.0

# Calculate SFS
vcftools --vcf /geneticData/variants/filtered_AF/combined/merged237/merged_combined.vcf --freq --out sfs_results

# Split last two columns to properly read into R
awk 'BEGIN {OFS="\t"} {split($5, a, ":"); print $1, $2, $3, $4, a[1], a[2]}' /geneticData/variants/filtered_AF/combined/merged237/sfs/sfs_results.frq > /geneticData/variants/filtered_AF/combined/merged237/sfs/sfs_preprocessed.frq