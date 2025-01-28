#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -t 96:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

CORES=8
BEDTOOLS="/home/dudekmf/local/bin/bedtools/bedtools"

module load bzip2

# Used to convert bams to fragment files
# From https://github.com/HYsxe/PRINT/issues/6

dir="../../../raw_data/brandon_liver_ATAC/namesorted_bams"
# bam_name="3782-BW-100"
# bam_file="$dir/$bam_name.namesorted.bam"

# bam file must be name sorted
for bam_file in $dir/*.namesorted.bam; do

    arrIN=(${bam_file//\// })
    arrIN=(${arrIN[-1]//.namesorted.bam/ })
    bam_name=${arrIN[0]}

    frags_name="frags/$bam_name.tsv"
    
    # if [ -f "$frags_name" ]; then
    #     echo "Skipping $bam_name because $frags_name already exists..."
    #     continue 
    # fi

    echo "Converting $bam_name to fragments file..."
    # Discard stderr to ignore warnings of unpaired reads
    $BEDTOOLS bamtobed -i $bam_file -bedpe 2> /dev/null | \
    awk -v OFS="\t" -v BARCODE="$bam_name" '{if($1==$4 && ($9=="+" || ($2==$5 && $3==$6))){print $1,$2+4,$6-5,BARCODE}}' | \
    sort --parallel=$CORES -S 40G  -k4,4 -k1,1 -k2,2n -k3,3n | \
    uniq -c | \
    awk -v OFS="\t" '{if($4>$3){print "chr"$2, $3, $4, $5, $1}}' > \
    $frags_name

    # Link file
    # mkdir "data/$bam_name"
    frags_path=$(realpath $frags_name)
    ln -sf $frags_path data/$bam_name/fragments.tsv

done
echo "Done"


