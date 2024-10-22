#!/bin/bash

#don't change this line
#SBATCH --account=def-koborlab

#Set max walltime for the job (hours:minutes:seconds)
#SBATCH --time=168:00:00

#Set memory requirements (don't request more than you need)
#SBATCH --cpus-per-task=1
#SBATCH --mem 100G

#Job name
#SBATCH --job-name=blat

#Output file
#SBATCH --output=blat.out

#Email notifications for: begin, end, fail, requeue, all
#SBATCH --mail-user=zdong@bcchr.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set variables

#your working directory
WORKING_DIR=/home/zdong01/scratch/zdong01/MONKEY/Conserved_Sequences_Check/Human

#location of your conda environments folder
#CONDA_PREFIX=$HOME/projects/def-koborlab/koborlab/sam_test/conda_envs

#location of personal/local R libraries in your conda environment
#export R_LIBS_USER=$HOME/projects/def-koborlab/sam_test/conda_envs/R3.5/lib/R/library/

#source the Kobor HPC shared software
source $HOME/projects/def-koborlab/koborlab/hpcenv_cedar.sh

#source your conda environment for personal packages
#source activate $HOME/projects/def-koborlab/koborlab/sam_test/conda_envs/R3.5

#load dependencies required to run R on compute canada
#module load nixpkgs/16.09
#module load gcc/7.3.0

#load R version 3.5.1
#R script
#module load r/3.5.1

#run your R script
bash blat.sh
#Rscript $HOME/projects/def-koborlab/koborlab/sam_test/ggplottest.R


