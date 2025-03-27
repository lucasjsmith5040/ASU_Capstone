#!/bin/bash
#SBATCH --job-name=fastqcAnalysis
#SBATCH --output=logs/fastqc_%j.out
#SBATCH --error=logs/fastqc_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

#Load fastqc module
module load fastqc-0.12.1-gcc-11.2.0

#Navigate to working directory
cd /scratch/dbihnam/lsc585/turtleProject/rawData

#Create directory for reports
mkdir -p FastqcOutput

#Run fastqc on all fastq.gz_1/2 files
fastqc -t 16 -o FastqcOutput *_1.fastq.gz *_2.fastq.gz

