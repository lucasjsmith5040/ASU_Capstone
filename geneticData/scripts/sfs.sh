#!/bin/bash

#Load vcftools
module load vcftools-0.1.14-gcc-11.2.0

#Calculate SFS
vcftools --vcf /scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/merged_combined.vcf --freq --out sfs_results

#Split last two columns to read into R
awk 'BEGIN {OFS="\t"} {split($5, a, ":"); print $1, $2, $3, $4, a[1], a[2]}' /scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/sfs/sfs_results.frq > /scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/sfs/sfs_preprocessed.frq