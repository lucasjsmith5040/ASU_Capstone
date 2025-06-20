#!/bin/bash
#SBATCH --job-name=fillTags
#SBATCH --output=logs/fillTags/fillTags_%j.out
#SBATCH --error=logs/fillTags/fillTags_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

## This script is intended to add allele frequency (AF) annotation to the INFO
## column of the vcf files

# Load bcftools
module load bcftools-1.14-gcc-11.2.0

# Set file paths
VCF_DIR="/geneticData/variants/unfiltered"
OUTPUT_DIR="/geneticData/variants/unfiltered_AF"
SAMPLES_TXT="/geneticData/rawData/Trimmed/paired/samples.txt"

# Create out directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Get the sample name using SLURM array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TXT")

# File paths
INPUT_VCF="${VCF_DIR}/${SAMPLE}.raw.vcf"
OUTPUT_VCF="${OUTPUT_DIR}/${SAMPLE}.raw_with_AF.vcf"

# Use fill-tags to add AF to the INFO field
bcftools +fill-tags "$INPUT_VCF" -o "$OUTPUT_VCF" -- -t AF
