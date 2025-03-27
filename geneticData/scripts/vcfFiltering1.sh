#!/bin/bash
#SBATCH --job-name=vcfFiltering1
#SBATCH --output=logs/vcfFiltering1/vcfFiltering_%j.out
#SBATCH --error=logs/vcfFiltering1/vcfFiltering_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

#Load the bcftools module
module load bcftools-1.14-gcc-11.2.0

#Define paths
VCF_DIR="/scratch/dbihnam/lsc585/turtleProject/variants/unfiltered"
FILTERED_DIR="/scratch/dbihnam/lsc585/turtleProject/variants/filteredTest"
SAMPLES_TXT="/scratch/dbihnam/lsc585/turtleProject/rawData/Trimmed/paired/samples.txt"

# Get sample name using SLURM array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TXT")

# File paths
INPUT_VCF="${VCF_DIR}/${SAMPLE}.raw.vcf"
OUTPUT_VCF="${FILTERED_DIR}/${SAMPLE}.filtered.vcf"

# Apply filters
bcftools filter -i 'QUAL > 30 && DP > 10 && F_MISSING < 0.1' \
  -o "$OUTPUT_VCF" "$INPUT_VCF"
