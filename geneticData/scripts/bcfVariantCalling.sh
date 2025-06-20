#!/bin/bash
#SBATCH --job-name=variantCalling
#SBATCH --output=logs/variantCalling/variantCalling_%j.out
#SBATCH --error=logs/variantCalling/variantCalling_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

## This script is intended to call variant SNPs and INDELs from individual BAM files
## This will result in individual VCF files per sample

# Load bcftools module
module load bcftools-1.14-gcc-11.2.0

# Define file paths
REF="/geneticData/references/ncbi_dataset/data/GCF_015237465.2/GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna"
BAM_DIR="/geneticData/alignedBam"
VCF_DIR="/geneticData/variants/unfiltered"
SAMPLES_TXT="/geneticData/rawData/Trimmed/paired/samples.txt"

# Get sample name from samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TXT")

# Generate the pileup and call variants
# Exclude variants with a mapping and base quality below 20
# Call SNPs and INDELs
bcftools mpileup -f "$REF" -q 20 -Q 20 -Ou "$BAM_DIR/$SAMPLE.sorted.bam" | \
bcftools call -mv -Ov -o "${VCF_DIR}/${SAMPLE}.raw.vcf"

