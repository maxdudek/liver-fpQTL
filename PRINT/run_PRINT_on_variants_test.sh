#!/bin/bash
#SBATCH --job-name="PRINT_variants"
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R
module load CUDA
module load hdf5
module load GCCcore/12.2.0

# source /mnt/isilon/sfgi/programs/miniconda3/bin/activate PRINT
# source ~/.bashrc
# conda activate PRINT
echo $PYTHONPATH

# python verify_h5py.py

Rscript run_PRINT_on_variants.R 48
