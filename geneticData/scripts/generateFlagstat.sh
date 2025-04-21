#!/bin/bash
#SBATCH --job-name=generateFlagstat
#SBATCH --output=logs/flagstat/flagstat_%j.out
#SBATCH --error=logs/flagstat/flagstat_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --array=1-268
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

#This script is intended to generate flagstat reports of all aligned BAM files
#in a provided directory
#Flagstat reports give alignment quality metrics

#Load Samtools module
module load samtools-1.16-gcc-11.2.0

#Navigate to BAM directory
cd /scratch/dbihnam/lsc585/turtleProject/alignedBam/GC_test

#Get sample names
samples=($(ls *.bam | grep -v '.bam.bai' | sed 's/.bam//'))
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

#Generate flagstat report
samtools flagstat ${sample}.bam > flagstatReports/${sample}_flagstat.txt