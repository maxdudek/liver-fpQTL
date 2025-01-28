#!/bin/bash
#SBATCH --job-name="splitSNP_cloglog"
#SBATCH -c 1
#SBATCH --mem=8GB
#SBATCH -t 24:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

source ~/.bashrc

FP_METHOD="PRINT_beta_cloglog"
HET_SAMPLES_FILE="het_samples/${FP_METHOD}_fpQTLs_het_samples.txt"

BAM_DIR="/mnt/isilon/sfgi/dudekm/raw_data/brandon_liver_ATAC/chr_renamed_bams"

# Create summary file
echo -e "args.output_prefix\tchrom\tpos\tref\talt\tref_count\talt_count" > allelic_fragments/${FP_METHOD}_splitSNP_summary.txt

while IFS=$'\t' read -r -a line
do
    rs=${line[0]}
    chrom=${line[1]}
    pos=${line[2]}
    ref=${line[3]}
    alt=${line[4]}
    het_samples=${line[6]}

    # Skip first line
    if [ "$rs" == "variant_id" ]; then
        continue
    fi

    echo $rs

    out_dir="split_bams/$FP_METHOD/${rs}_$chrom:$pos"
    mkdir -p $out_dir
    echo "$out_dir" > $out_dir/log.txt

    # Iterate over heterozygous samples
    while IFS=',' read -ra ADDR; do
    for sample in "${ADDR[@]}"; do
        bam="$BAM_DIR/$sample.bam"

        echo "" >> $out_dir/log.txt
        python splitSNP.py $bam $out_dir/${sample}_$rs $chrom:$pos:$ref:$alt &>> $out_dir/log.txt
    done
    done <<< "$het_samples"

    # break
done < $HET_SAMPLES_FILE



