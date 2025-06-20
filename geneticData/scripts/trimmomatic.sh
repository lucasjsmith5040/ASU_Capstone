#!/bin/bash
#SBATCH --job-name=trimmomatic_analysis
#SBATCH --output=logs/trimmomatic/trimmomatic_%j.out
#SBATCH --error=logs/trimmomatic/trimmomatic_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

## This script is intended to trim raw FASTQ reads, removing adapters and low quality
## sequences to prepare for downstream genomic alignment

# Load in Trimmomatic module
module load trimmomatic-0.39-gcc-12.1.0

# Navigate to directory
cd /geneticData/rawData

# Make sure trimmed directory exists
mkdir -p Trimmed

# Obtain sample names
samples=($(ls *_1.fastq.gz | sed 's/_1.fastq.gz//'))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

# Run Trimmomatic
# Remove Illumina adapters using TruSeq3-PE.fa with max 2 mismatches, palindrome threshold 30, simple clip threshold 10
# Trim low-quality bases from the start of reads (below Q3)
# Trim low-quality bases from the end of reads (below Q3)
# Trim when average quality in a 4-base sliding window falls below Q20
# Discard reads shorter than 50 bases after trimming
trimmomatic PE -threads 4 \
  ${sample}_1.fastq.gz ${sample}_2.fastq.gz \
  Trimmed/${sample}_1_paired.fastq.gz Trimmed/${sample}_1_unpaired.fastq.gz \
  Trimmed/${sample}_2_paired.fastq.gz Trimmed/${sample}_2_unpaired.fastq.gz \
  ILLUMINACLIP:/scratch/dbihnam/lsc585/turtleProject/references/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
  