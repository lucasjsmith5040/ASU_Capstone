#!/bin/bash
#SBATCH --job-name=vcfFiltering_AF
#SBATCH --output=logs/vcfFiltering1/vcfFilteringAF_%j.out
#SBATCH --error=logs/vcfFiltering1/vcfFilteringAF_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

#This script was used to filter the VCF files created from the previous
#variant calling step, as well as separate out heterozygous and homozygous variants
#into their own separate vcf files for possible future analyses
#This set of filtering conditions is more stringent than 'vcfFiltering1.sh',
#and was used for the final analysis of the data

# Load the bcftools module
module load bcftools-1.14-gcc-11.2.0

#Define paths
VCF_DIR="/scratch/dbihnam/lsc585/turtleProject/variants/unfiltered_AF"
FILTERED_DIR="/scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF"
SAMPLES_TXT="/scratch/dbihnam/lsc585/turtleProject/rawData/Trimmed/paired/samples.txt"

#Get sample name using SLURM array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TXT")

#File paths
INPUT_VCF="${VCF_DIR}/${SAMPLE}.raw_with_AF.vcf"
OUTPUT_VCF="${FILTERED_DIR}/${SAMPLE}.filtered.vcf"

#Apply vcf filters
#phred QUAL > 30 (99.9% confidence)
#sequencing depth > 30 
#missing genotypes < 10%
bcftools filter -i 'QUAL > 30 && DP > 30 && F_MISSING < 0.1' \
  -o "$OUTPUT_VCF" "$INPUT_VCF"

#Create new copies of the filtered VCF file from above
#One file with only homozygous variants, another with heterozygous
bcftools view -i 'GT="1/1" || GT="0/0"' -o "${FILTERED_DIR}/${SAMPLE}.homozygous.vcf" "$OUTPUT_VCF"
bcftools view -i 'GT="0/1" || GT="1/0"' -o "${FILTERED_DIR}/${SAMPLE}.heterozygous.vcf" "$OUTPUT_VCF"

#Verify all samples are filtered
echo "Filtering and splitting completed for sample: $SAMPLE"
