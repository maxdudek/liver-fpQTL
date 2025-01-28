#!/bin/bash
#SBATCH --job-name="get_allelic_fragments_cloglog"
#SBATCH -c 1
#SBATCH --mem=8GB
#SBATCH -t 24:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

source ~/.bashrc
module load bzip2

BEDTOOLS="/home/dudekmf/local/bin/bedtools/bedtools"

# Create fragments file
FP_METHOD="PRINT_beta_cloglog"
HET_SAMPLES_FILE="het_samples/${FP_METHOD}_fpQTLs_het_samples.txt"
OUTFILE="allelic_fragments/${FP_METHOD}_allelic_fragments.txt"
echo -e "rsID\tsample\tallele\tinsertion_5prime\tinsertion_3prime\tcount" > $OUTFILE

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

    # Iterate over heterozygous samples
    while IFS=',' read -ra ADDR; do
    for sample in "${ADDR[@]}"; do
        ref_bam="$out_dir/${sample}_$rs.ref.$ref.bam"
        alt_bam="$out_dir/${sample}_$rs.alt.$alt.bam"

        # Namesort bams
        samtools sort -n -o $ref_bam.namesorted $ref_bam
        samtools sort -n -o $alt_bam.namesorted $alt_bam

        # Write ref allele fragments
        $BEDTOOLS bamtobed -i $ref_bam.namesorted -bedpe 2> /dev/null | \
        awk -v OFS="\t" -v RSID="$rs" -v SAMPLE="$sample" '{if($1==$4 && ($9=="+" || ($2==$5 && $3==$6))){print $2+4,$6-5,RSID,SAMPLE,"ref"}}' | \
        sort -k1,1n -k2,2n | \
        uniq -c | \
        awk -v OFS="\t" '{if($3>$2){print $4, $5, $6, $2, $3, $1}}' >> \
        $OUTFILE

        # Write alt allele fragments
        $BEDTOOLS bamtobed -i $alt_bam.namesorted -bedpe 2> /dev/null | \
        awk -v OFS="\t" -v RSID="$rs" -v SAMPLE="$sample" '{if($1==$4 && ($9=="+" || ($2==$5 && $3==$6))){print $2+4,$6-5,RSID,SAMPLE,"alt"}}' | \
        sort -k1,1n -k2,2n | \
        uniq -c | \
        awk -v OFS="\t" '{if($3>$2){print $4, $5, $6, $2, $3, $1}}' >> \
        $OUTFILE
        
    done
    done <<< "$het_samples"

    # break
done < $HET_SAMPLES_FILE
