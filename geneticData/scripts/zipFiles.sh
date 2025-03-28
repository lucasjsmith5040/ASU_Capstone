#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=compress_vcf
#SBATCH --output=logs/zipFiles/compress_vcf_%j.out
#SBATCH --error=logs/zipFiles/compress_vcf_%j.err
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

#Load pigz
module load pigz-2.6-gcc-11.2.0

#Load bgzip
module load htslib-1.14-gcc-11.2.0

#Check for input arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path/to/you/files> <fileExtension_NO_DOT>"
    exit 1
fi

#Assign file location and extension variables
DIR=$1
EXT=$2

#Navigate to user-specified directory
cd "$DIR" || exit 1

#Zip files checking for ext type
if [ "$EXT" == "vcf" ]; then
  #Fix my mistake of zipping the vcfs with pigz first lol
  for file in *.$EXT.gz; do
      if [ -f "$file" ]; then
        pigz -d "$file"
      fi
  done
  #Zip vcf with bgzip
  for file in *.$EXT; do
    bgzip -f "$file"
  done
else
  pigz -p 16 *.$EXT
fi
