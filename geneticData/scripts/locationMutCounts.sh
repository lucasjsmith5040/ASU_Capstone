#!/bin/bash

## Results from this script were not used in the final paper
## Sample size and average number of mutations per sample were plotted
## on a map of the region with ArcGIS, but the map was left out of the final paper

# Navigate to directory
cd /geneticData/variants/filtered_AF/combined/merged237/location

# Match all vcfs
vcf_files=(merged_*_samples.vcf.gz)

# Set output file
output_file="sampleLocationMutCounts.tsv"

# Write output header
echo -e "Region\tAverage_Variants_Per_Sample\tNum_Samples" > "$output_file"

# Loop through all vcfs and get mutation count statistics
for vcf in "${vcf_files[@]}"; do
  region=$(basename "$vcf" | sed 's/\.[^.]*$//')
  num_samples=$(bcftools query -l "$vcf" | wc -l)
  avg=$(bcftools query -f '[%GT\t]\n' "$vcf" | \
  awk '
  {
    for(i=1; i<=NF; i++) {
      if($i != "./." && $i != "." && $i != ".|.") count[i]++
    }
  }
  END {
    total=0;
    for(i in count) {
      total += count[i];
    }
    if(length(count) > 0)
      print total/length(count);
    else
      print 0;
  }')

# Append results to output
  echo -e "$region\t$avg\t$num_samples" >> "$output_file"
done
