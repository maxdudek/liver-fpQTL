#!/bin/bash
#SBATCH --job-name="extract_cutsites"
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.2.3

Rscript extract_cutsites_chr17.R
