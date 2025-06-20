#!/bin/bash

## This script is intended to merge the vcfs from a user-specified subsetted list
## based on the sample source location
## Used for site-specific SFS

# Check if the sample file argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <sample_file.txt>"
    exit 1
fi

# Use user-specified location file
sample_file=$1

# Initialize the list of VCFs to be merged
vcf_list=""

# Loop through the sample IDs in the provided sample file
while read sample; do
    vcf_list="$vcf_list ${sample}.filtered.vcf.gz"
done < "$sample_file"

# Merge vcfs
bcftools merge $vcf_list -o merged237/location/merged_${sample_file%.txt}.vcf.gz -Oz
