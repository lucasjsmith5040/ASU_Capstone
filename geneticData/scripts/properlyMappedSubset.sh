#!/bin/bash

## This script is intended to parse the flagstat reports of each sample, and generate
## a new list of samples that have a properly mapped reads percentage >70%
## This yielded a list of 237 samples, down from 268

# Define the input file and threshold
INPUT_FILE="/geneticData/alignedBam/flagstatReports/flagstat_summary.csv"
THRESHOLD=70

# Output file with valid samples
VALID_SAMPLES_FILE="valid_samples.txt"

# Process the CSV file and calculate the percentage of properly paired reads
awk -F ',' 'NR > 1 { 
    properly_paired_reads = $3 
    total_reads = $2 
    percentage = (properly_paired_reads / total_reads) * 100
    if (percentage >= 70) print $1
}' "$INPUT_FILE" > "$VALID_SAMPLES_FILE"

echo "Valid samples with > 70% properly paired reads are saved in $VALID_SAMPLES_FILE"
