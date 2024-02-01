#!/bin/bash

#SBATCH --qos=normal
#SBATCH --job-name MIM2_8450
#SBATCH --error /lustre/project/svanbael/bolivar/Mimulus_sequences/MIM2_8450.error          
#SBATCH --output /lustre/project/svanbael/bolivar/Mimulus_sequences/MIM2_8450.output  
#SBATCH --time=23:00:00
#SBATCH --mem=256000 #Up to 256000
#SBATCH --nodes=4               # One is enough. When running MPIs,anywhere from 2-4 should be good.
#SBATCH --ntasks-per-node=1    # Number of Tasks per Node
#SBATCH --cpus-per-task=20      # Number of threads per task (OMP threads)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=baponterolon@tulane.edu
#SBATCH --partition=centos7    #This is important to run the latest software versions

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK #$SLURM_CPUS_PER_TASK stores whatever value we assign to cpus-per-task, and is therefore our candidate for passing to OMP_NUM_THREADS

##Modules/Singularity
module load anaconda3/2020.07
source activate virtual_env

#Run R script
R CMD BATCH dada2_mim_pipeline.R

module purge