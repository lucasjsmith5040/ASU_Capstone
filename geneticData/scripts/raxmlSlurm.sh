#!/bin/bash
#SBATCH --job-name=raxml_phylogeneticTree_MPI
#SBATCH --output=logs/raxml/raxml_%j.out
#SBATCH --error=logs/raxml/raxml_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --mail-user=dbihnam@asu.edu
#SBATCH --mail-type=ALL

#Load MPI-enabled raxml-ng
module load raxml-ng-1.1.0-gcc-11.2.0 

#Run raxml-ng-mpi for distributed processing
raxml-ng-mpi --msa /scratch/dbihnam/lsc585/turtleProject/variants/filtered_AF/combined/merged237/phylogenetics/merged_combined.min4.phy --model GTR+G --prefix output_tree --seed 12345
