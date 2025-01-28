#!/bin/bash
#SBATCH --job-name="extract_cutsites"
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH -t 3-00:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.4.0

for FP_METHOD in PRINT_beta_cloglog; do
    mkdir -p cutsites/$FP_METHOD
    echo "Extracting cutsites for $FP_METHOD..."
    echo ""
    Rscript extract_cutsites.R $FP_METHOD
    echo ""
done
