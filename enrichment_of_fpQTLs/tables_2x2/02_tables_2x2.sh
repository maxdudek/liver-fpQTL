#!/bin/bash
#SBATCH --job-name="tables_2x2"
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.4.0

Rscript tables_2x2.R

