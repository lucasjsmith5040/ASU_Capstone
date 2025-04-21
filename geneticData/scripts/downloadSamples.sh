#!/bin/bash

#This is a short script to download all the raw NCBI samples with sra ids in
#the provided text file and convert them to FASTQ format

#Set output directory
OUTPUT="/scratch/dbihnam/lsc585"

#Read each line of the input file
while IFS= read -r sra_id; do
echo "Downloading and converting $sra_id..."
#Convert SRA file to FASTQ format
fasterq-dump $sra_id --outdir $OUTPUT
#Provide input file
done < sra_ids1.txt

