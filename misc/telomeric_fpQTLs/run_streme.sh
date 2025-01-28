#!/bin/bash
#SBATCH --job-name="streme"
#SBATCH -c 1
#SBATCH --mem=32GB
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

source ~/.bashrc

OUT_DIR="streme_output"

# STREME parameters
MINW=6
MAXW=20

HOURS=6
SECS=$(($HOURS*3600))

echo "SECS = $SECS"
OPTS="--minw $MINW --maxw $MAXW --verbosity 2 --time $SECS"

# Telomeric fpQTLs vs. non-telomeric fpQTLs
# echo "Telomeric fpQTLs vs. non-telomeric fpQTLs"
# streme $OPTS \
#     --oc ${OUT_DIR}/telomeric_fpQTLs_vs_non_telomeric_fpQTLs \
#     --p seqs/telomeric_fpQTLs.fasta \
#     --n seqs/non_telomeric_fpQTLs.fasta

# Telomeric fpQTLs vs other variants
echo "Telomeric fpQTLs vs other variants"
streme $OPTS \
    --oc ${OUT_DIR}/telomeric_fpQTLs_vs_other_variants \
    --p seqs/telomeric_fpQTLs.fasta \
    --n seqs/non_telomericfpQTL_variants.fasta

# fpQTLs vs non_fpQTLs
echo "fpQTLs vs non_fpQTLs"
streme $OPTS \
    --oc ${OUT_DIR}/fpQTLs_vs_non_fpQTLs \
    --p seqs/all_fpQTLs.fasta \
    --n seqs/non_fpQTLs.fasta

# fpQTLs vs none
echo "fpQTLs vs none"
streme $OPTS \
    --oc ${OUT_DIR}/fpQTLs_vs_none \
    --p seqs/all_fpQTLs.fasta

