#!/bin/bash
#SBATCH --job-name="subset_bams"
#SBATCH -c 1
#SBATCH --mem=32GB
#SBATCH -t 48:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

source ~/.bashrc
module load bzip2

BAM_DIR="/mnt/isilon/sfgi/dudekm/raw_data/brandon_liver_ATAC/chr_renamed_bams"

for bam_file in $BAM_DIR/*.bam; do

    arrIN=(${bam_file//\// })
    arrIN=(${arrIN[-1]//.bam/ })
    bam_name=${arrIN[0]}

    outfile="bams/${bam_name}_chr17_inversion.bam"

    echo $bam_name 
    # samtools view -o $outfile $bam_file chr17:45000000-47000000
    samtools index $outfile

done

