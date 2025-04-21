#!/bin/bash

#This script is used to extract a user-specified number of principal components (PCs)
#from a merged vcf file

#Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <Number_of_PCs>"
    exit 1
fi

#Set the number of PCs from the argument (Must be a positive integer)
NPCS=$1

#Load plink module
module load plink

#Create $NPCS directory
mkdir -p /scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/pca/$NPCS

#Convert merged VCF to genotype format
plink \
  --vcf /scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/merged_combined.vcf \
  --make-bed \
  --out /scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/pca/$NPCS/merged_combined \
  --allow-extra-chr

#Run PCA using the specified number of principal components
plink \
  --bfile /scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/pca/$NPCS/merged_combined \
  --pca "$NPCS" \
  --out /scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/pca/$NPCS/pcaResults \
  --allow-extra-chr
