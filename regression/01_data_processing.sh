#!/bin/bash
#SBATCH --job-name="data_processing"
#SBATCH -c 1
#SBATCH --mem=200G
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.4.0

# FP_METHOD="TOBIAS"
# FP_METHOD="TOBIAS_cpm_normalized"

for FP_METHOD in PRINT_no_gaussian; do
    echo "Processing data for $FP_METHOD..."
    echo ""
    Rscript data_processing.R $FP_METHOD
done




