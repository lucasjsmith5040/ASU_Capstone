#!/bin/bash
#SBATCH --job-name=jointVariantCalling
#SBATCH --output=logs/variantCalling/jointVariantCalling_%j.out
#SBATCH --error=logs/variantCalling/jointVariantCalling_%j.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

# Load bcftools module
module load bcftools-1.14-gcc-11.2.0

# Define file paths
REF="/scratch/dbihnam/lsc585/turtleProject/references/ncbi_dataset/data/GCF_015237465.2/GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna"
BAM_DIR="/scratch/dbihnam/lsc585/turtleProject/alignedBam"
VCF_DIR="/scratch/dbihnam/lsc585/turtleProject/variants/population/unfiltered"
VALID_SAMPLES_FILE="/scratch/dbihnam/lsc585/turtleProject/alignedBam/flagstatReports/valid_samples.txt"

# Create a list of BAM file paths from valid_samples.txt
BAM_LIST=$(cat "$VALID_SAMPLES_FILE" | while read SAMPLE; do echo "$BAM_DIR/$SAMPLE"; done)

# Generate pileup and call variants for all valid samples together
bcftools mpileup -f "$REF" -q 20 -Q 20 -Ou $BAM_LIST | \
bcftools call -mv -Ov -o "${VCF_DIR}/population_raw.vcf"

# Optional: Filter VCF based on quality, depth, etc.
bcftools filter -i 'QUAL>30 && DP>10' "${VCF_DIR}/population_raw.vcf" > "${VCF_DIR}/population_filtered.vcf"
