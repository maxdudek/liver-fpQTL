#!/bin/bash
#SBATCH --job-name="get_motif_strength"
#SBATCH -c 1
#SBATCH --mem=500G
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.2.3

Rscript get_motif_strength.R 5e-5
Rscript get_motif_strength.R 5e-4

