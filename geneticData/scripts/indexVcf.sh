#!/bin/bash

## This script is intended to index all the vcf files in a user-specified directory

# Store user-specified directory
DIR=$1

# Load bcftools
module load bcftools-1.14-gcc-11.2.0

# Check for directory argument
if [ -z "$1" ]; then
  echo "Usage: $0 <path/to/vcf/files>"
  exit 1
fi

# Loop through directory and index vcf
for vcf in "$DIR"/*.vcf.gz; do
  bcftools index -t "$vcf"
done
