#!/bin/bash
#SBATCH --job-name="summarize_motif_results"
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.4.0

Rscript summarize_motif_results.R p=5e-4
Rscript summarize_motif_results.R p=5e-5

