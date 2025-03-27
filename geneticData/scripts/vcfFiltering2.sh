#!/bin/bash
#SBATCH --job-name=vcfFiltering2
#SBATCH --output=logs/vcfFiltering/vcfFiltering_%j.out
#SBATCH --error=logs/vcfFiltering/vcfFiltering_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

# Load the bcftools module
module load bcftools-1.14-gcc-11.2.0

# Define paths
VCF_DIR="/scratch/dbihnam/lsc585/turtleProject/variants/unfiltered"
FILTERED_DIR="/scratch/dbihnam/lsc585/turtleProject/variants/filtered2"
SAMPLES_TXT="/scratch/dbihnam/lsc585/turtleProject/rawData/Trimmed/paired/samples.txt"

# Get sample name using SLURM array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TXT")

# File paths
INPUT_VCF="${VCF_DIR}/${SAMPLE}.raw.vcf"
OUTPUT_VCF="${FILTERED_DIR}/${SAMPLE}.filtered.vcf"

# Copy the raw VCF file to the filtered directory to avoid modifying the original
cp "$INPUT_VCF" "$OUTPUT_VCF"

# Calculate allele frequencies and add AF tag to INFO field
bcftools +fill-tags "$OUTPUT_VCF" -o "$OUTPUT_VCF"

# Apply filters: QUAL > 30, DP > 10, AF between 0.05 and 0.50
bcftools filter -i 'QUAL > 30 && DP > 10 && AF > 0.05 && AF < 0.50' \
  -o "$OUTPUT_VCF" "$OUTPUT_VCF"

# Zip and index the filtered VCF file
bgzip "$OUTPUT_VCF"
tabix -p vcf "${OUTPUT_VCF}.gz"
