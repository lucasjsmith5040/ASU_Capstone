#!/bin/bash
#SBATCH --job-name=fillTags
#SBATCH --output=logs/fillTags/fillTags_%j.out
#SBATCH --error=logs/fillTags/fillTags_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

# Load the necessary bcftools module
module load bcftools-1.14-gcc-11.2.0

# Define paths
VCF_DIR="/scratch/dbihnam/lsc585/turtleProject/variants/unfiltered"
OUTPUT_DIR="/scratch/dbihnam/lsc585/turtleProject/variants/unfiltered_AF"
SAMPLES_TXT="/scratch/dbihnam/lsc585/turtleProject/rawData/Trimmed/paired/samples.txt"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Get the sample name using SLURM array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TXT")

# File paths
INPUT_VCF="${VCF_DIR}/${SAMPLE}.raw.vcf"
OUTPUT_VCF="${OUTPUT_DIR}/${SAMPLE}.raw_with_AF.vcf"

# Apply the fill-tags to add AF to the INFO field
bcftools +fill-tags "$INPUT_VCF" -o "$OUTPUT_VCF" -- -t AF
