#!/bin/bash
#SBATCH --job-name="simulations"
#SBATCH -c 32
#SBATCH --mem=400G
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.4.0

Rscript simulations.R 32
