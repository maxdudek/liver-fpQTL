#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

export HDF5_DISABLE_VERSION_CHECK=2

module load R/4.4.0
module load CUDA
module load hdf5
module load GCCcore/12.2.0

Rscript run_PRINT.R PRINT_beta_cloglog
