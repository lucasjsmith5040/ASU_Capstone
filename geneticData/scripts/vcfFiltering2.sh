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
bcftools filter -i 'QUAL > 30 && DP > 30 && F_MISSING < 0.1' \
  -o "$OUTPUT_VCF" "$INPUT_VCF"

#Split into homozygous and heterozygous VCF files
bcftools view -i 'GT="1/1" || GT="0/0"' -o "${FILTERED_DIR}/${SAMPLE}.homozygous.vcf" "$OUTPUT_VCF"
bcftools view -i 'GT="0/1" || GT="1/0"' -o "${FILTERED_DIR}/${SAMPLE}.heterozygous.vcf" "$OUTPUT_VCF"

echo "Filtering and splitting completed for sample: $SAMPLE"
