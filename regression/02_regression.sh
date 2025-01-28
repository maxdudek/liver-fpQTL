#!/bin/bash
#SBATCH --job-name="regression_32"
#SBATCH -c 32
#SBATCH --mem=800G
#SBATCH -t 10-0:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.4.0

# Rscript beta_regression.R 32

Rscript linear_regression.R PRINT_no_gaussian 32

# for FP_METHOD in PRINT_no_gaussian; do
#     echo ""
#     Rscript linear_regression.R $FP_METHOD 16
#     echo ""
# done

# Get residuals
# for FP_METHOD in PRINT_no_gaussian; do
#     echo ""
#     Rscript get_residuals.R $FP_METHOD 16
#     echo ""
# done

