#!/bin/bash
#SBATCH --job-name="plot_cutsites"
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 72:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.4.0

for FP_METHOD in PRINT_beta_cloglog PRINT_beta_logit; do
    mkdir -p figures/$FP_METHOD
    echo "Plotting cutsites for $FP_METHOD..."
    echo ""
    Rscript plot_cutsites.R $FP_METHOD
    echo ""
done
