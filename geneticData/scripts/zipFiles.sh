#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=compress_fastq
#SBATCH --output=compress_fastq_%j.out
#SBATCH --error=compress_fastq_%j.err
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

#Load pigz
module load pigz-2.6-gcc-11.2.0

#Navigate to directory
cd /scratch/dbihnam/lsc585/turtleProject/rawData

#Compress all FASTQ files
pigz -p 16 *.fastq
