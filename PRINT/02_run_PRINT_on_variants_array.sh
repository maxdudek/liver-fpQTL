#!/bin/bash
#SBATCH --job-name="PRINT_variants"
#SBATCH -c 2
#SBATCH --mem=32G
#SBATCH -t 96:00:00
#SBATCH -o job_out/array3/%x.%A_%a.out 
#SBATCH -e job_out/array3/%x.%A_%a.out
#SBATCH --array=1-189

# if [[ "$SLURM_ARRAY_TASK_ID" == "48" ]]; then
#     echo "DUDEK: skipping 48"
#     exit
# fi

module load R
module load CUDA
module load hdf5
module load GCCcore/12.2.0

Rscript run_PRINT_on_variants.R $SLURM_ARRAY_TASK_ID
