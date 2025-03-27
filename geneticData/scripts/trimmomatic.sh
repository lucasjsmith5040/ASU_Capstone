#!/bin/bash
#SBATCH --job-name=trimmomatic_analysis
#SBATCH --output=logs/trimmomatic/trimmomatic_%j.out
#SBATCH --error=logs/trimmomatic/trimmomatic_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

#Load in Trimmomatic module
module load trimmomatic-0.39-gcc-12.1.0

#Navigate to directory
cd /scratch/dbihnam/lsc585/turtleProject/rawData

#Make sure trimmed directory exists
mkdir -p Trimmed

#Obtain sample names
samples=($(ls *_1.fastq.gz | sed 's/_1.fastq.gz//'))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

#Run Trimmomatic
trimmomatic PE -threads 4 \
  ${sample}_1.fastq.gz ${sample}_2.fastq.gz \
  Trimmed/${sample}_1_paired.fastq.gz Trimmed/${sample}_1_unpaired.fastq.gz \
  Trimmed/${sample}_2_paired.fastq.gz Trimmed/${sample}_2_unpaired.fastq.gz \
  ILLUMINACLIP:/scratch/dbihnam/lsc585/turtleProject/references/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
  