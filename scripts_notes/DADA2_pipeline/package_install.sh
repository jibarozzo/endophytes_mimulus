#!/bin/bash

#SBATCH --qos=normal
#SBATCH --job-name package_install
#SBATCH --error /lustre/project/svanbael/bolivar/package_install.error          
#SBATCH --output /lustre/project/svanbael/bolivar/package_install.output  
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --partition=centos7

echo "Start job"

#Modules/Singularity
module load anaconda3/2020.07
source activate virtual_env

#Start R
R CMD BATCH dada2_package_install.R 

module purge

echo "end"


