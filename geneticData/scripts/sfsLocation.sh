#!/bin/bash

## This script is intended to calculate the site frequency spectrum (SFS) for 
## site-specific merged VCF files

# Load vcftools
module load vcftools-0.1.14-gcc-11.2.0

# Directory containing VCF files
vcf_dir="/geneticData/variants/filtered_AF/combined/merged237/location"

# Loop through each VCF file in the directory (assuming they are .vcf.gz files)
for vcf_file in $vcf_dir/*.vcf.gz; do
    base_filename=$(basename "$vcf_file" .vcf.gz)
    # Calculate SFS for each VCF file
    vcftools --gzvcf "$vcf_file" --freq --out "$vcf_dir/sfs_results_$base_filename"
    # Split the last two columns for R, for each VCF
    awk 'BEGIN {OFS="\t"} {split($5, a, ":"); print $1, $2, $3, $4, a[1], a[2]}' "$vcf_dir/sfs_results_$base_filename.frq" > "$vcf_dir/../sfs/location/sfs_preprocessed_$base_filename.frq"
done
