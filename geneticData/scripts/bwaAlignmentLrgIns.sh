#!/bin/bash
#SBATCH --job-name=bwaAlignmentLrgIns
#SBATCH --output=logs/bwaAlignment/bwaAlignment_%j.out
#SBATCH --error=logs/bwaAlignment/bwaAlignment_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

##This script was tried for aligning with an assumption of a large sequencing insert
##Generated bams not used in the downstream analysis

#Load BWA alignment tool
module load bwa-0.7.17-gcc-12.1.0

#Load Samtools module
module load samtools-1.16-gcc-11.2.0

#Define file paths
REF="/scratch/dbihnam/lsc585/turtleProject/references/ncbi_dataset/data/GCF_015237465.2/GCF_015237465.2_rCheMyd1.pri.v2_genomic.fna"
FASTQ_DIR="/scratch/dbihnam/lsc585/turtleProject/rawData/Trimmed/paired"
OUT_DIR="/scratch/dbihnam/lsc585/turtleProject/alignedBam/LrgIns"
SAMPLES_TXT="/scratch/dbihnam/lsc585/turtleProject/rawData/Trimmed/paired/samples.txt"

#Get sample name from samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TXT")

#Run alignment
#Optimized for larger insert sizes
#P keeps all read pairs to include ones with large fragment sizes
#a aligns full-length fragments instead of treating them as separate reads
#T 30 adjusts the alignment score threshold to improve handling of mismatches
#M marks shorter split hits as secondary
bwa mem -t 8 -P -a -T 30 -M "$REF" "${FASTQ_DIR}/${SAMPLE}_1_paired.fastq.gz" "${FASTQ_DIR}/${SAMPLE}_2_paired.fastq.gz" | \
  samtools view -bS | samtools sort -o "${OUT_DIR}/${SAMPLE}.sorted.bam"
  
#Index bam
samtools index "${OUT_DIR}/${SAMPLE}.sorted.bam"