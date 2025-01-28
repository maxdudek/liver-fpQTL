#!/bin/bash
#SBATCH --job-name="plotting_results"
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.4.0

# for FP_METHOD in PRINT_beta_logit  PRINT_beta_probit PRINT_beta_cloglog PRINT_beta_logit_phi PRINT_beta_probit_phi PRINT_beta_cloglog_phi; do
for FP_METHOD in PRINT_beta_cloglog PRINT_beta_logit; do
    echo "Plotting results for $FP_METHOD..."
    echo ""
    Rscript plotting_results.R $FP_METHOD
    echo ""
done
