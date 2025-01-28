#!/bin/bash

module load bzip2
module load R/4.4.0
PYTHONPATH=""

source ~/.bashrc
for i in $(seq ${CONDA_SHLVL}); do
    conda deactivate
done
conda activate /mnt/isilon/sfgi/programs/miniconda3/envs/ldsc

LDSC="/mnt/isilon/sfgi/programs/ldsc"

# Test LDSC environment
which python
# python $LDSC/make_annot.py -h
# python $LDSC/ldsc.py -h
# exit


# Khanh hg38
cdir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsc
BFILE=$cdir/plink_files_hg38/1000G.EUR.hg38
FREQ=$cdir/plink_files_hg38/1000G.EUR.hg38.
BASELINE=$cdir/baselineLD_v2.2_hg38/baselineLD.
WEIGHTS=$cdir/weights_hg38/weights.hm3_noMHC.
PRINT_SNPS=$cdir/list.txt

OCR_LDSCORE="../ldsr_discrete/annotations/brandon_liver_ocr/ld_score/brandon_liver_ocr."

BED_IN="annotations/$NAME/$NAME.bed"
LD_SCORE="annotations/$NAME/ld_score"

mkdir -p annotations/$NAME/ld_score
mkdir -p annotations/$NAME/out_baseline
mkdir -p annotations/$NAME/out_ocr
mkdir -p annotations/$NAME/out_ocr_baseline

echo "Running LDSR for functional annotation '$NAME'"

# Create annotation+LD files
if [ -f "$LD_SCORE/$NAME.22.l2.ldscore.gz" ]; then
    echo "Skipping annotation creation for $NAME because $LD_SCORE/$NAME.22.l2.ldscore.gz already exists..." 
else
    # Bed to Annot
    for j in {1..22}; do
    Rscript make_ldsc_continuous_annot.R \
        $BED_IN \
        $BFILE.$j.bim \
        $LD_SCORE/$NAME.$j.annot.gz "full-annot"
    done

    #Annot to LD files
    for j in {1..22}; do
    python $LDSC/ldsc.py \
        --l2 \
        --bfile $BFILE.$j \
        --ld-wind-cm 1 \
        --annot $LD_SCORE/$NAME.$j.annot.gz \
        --out $LD_SCORE/$NAME.$j \
        --print-snps $PRINT_SNPS
    done
fi

# Create quantile M files
BASELINE_QM="annotations/$NAME/out_baseline/baselineLD_$NAME.q5.M"
OCR_QM="annotations/$NAME/out_ocr/OCR_$NAME.q5.M"
OCR_BASELINE_QM="annotations/$NAME/out_ocr_baseline/OCR_baselineLD_$NAME.q5.M"
if [ -f "$OCR_BASELINE_QM" ]; then
    echo "Skipping QM creation for $NAME because $OCR_BASELINE_QM already exists..." 
else
    # perl $LDSC/ContinuousAnnotations/quantile_M.pl \
    #         --ref-annot-chr $LD_SCORE/$NAME.,$BASELINE \
    #         --frqfile-chr $FREQ \
    #         --annot-header "ANNOT" \
    #         --nb-quantile 5 \
    #         --maf 0.05 \
    #         --out $BASELINE_QM
            
    # perl $LDSC/ContinuousAnnotations/quantile_M.pl \
    #         --ref-annot-chr $OCR_LDSCORE,$LD_SCORE/$NAME. \
    #         --frqfile-chr $FREQ \
    #         --annot-header "ANNOT" \
    #         --nb-quantile 5 \
    #         --maf 0.05 \
    #         --out $OCR_QM
            
    perl $LDSC/ContinuousAnnotations/quantile_M.pl \
            --ref-annot-chr $OCR_LDSCORE,$BASELINE,$LD_SCORE/$NAME. \
            --frqfile-chr $FREQ \
            --annot-header "ANNOT" \
            --nb-quantile 5 \
            --maf 0.05 \
            --out $OCR_BASELINE_QM
fi



# Run LDSR, liver curated sumstats

# Option 1: selected traits
# SUMSTATS="/mnt/isilon/sfgi/dudekm/raw_data/GWAS/sumstats/munged_sumstats"
# for trait in ALT ALP GGT SCZ BIP lipids_TG lipids_TC lipids_HDL lipids_nonHDL lipids_LDL NAFLD BMI T2D; do
#     sumstat=$SUMSTATS/$trait.sumstats.gz

# Option 2: all downloaded traits
# SUMSTATS="/mnt/isilon/sfgi/dudekm/raw_data/GWAS/sumstats/munged_sumstats"
# for sumstat in "$SUMSTATS"/*.sumstats.gz; do
    # */NAFLD.sumstats.gz --> NAFLD
    # arrIN=(${sumstat//\// })
    # arrIN=(${arrIN[-1]//.sumstats./ })
    # trait=${arrIN[0]}

# Option 3: Price lab traits
SUMSTATS="/mnt/isilon/sfgi/dudekm/raw_data/GWAS/sumstats/alkesgroup_sumstats"
for sumstat in "$SUMSTATS"/*.sumstats.gz; do
    # */NAFLD.sumstats.gz --> NAFLD
    arrIN=(${sumstat//\// })
    arrIN=(${arrIN[-1]//.sumstats./ })
    trait=${arrIN[0]}


    echo "Running regression on $sumstat..."

    # Run annot vs. baseline
    # echo "annot vs. baseline"
    # OUT="annotations/$NAME/out_baseline/$trait.$NAME"
    # if [ -f "$OUT.q5.txt" ]; then
    #     echo "Skipping because $OUT.q5.txt already exists..."
    # else
    #     python $LDSC/ldsc.py \
    #         --h2 $sumstat \
    #         --ref-ld-chr $LD_SCORE/$NAME.,$BASELINE \
    #         --w-ld-chr $WEIGHTS \
    #         --overlap-annot \
    #         --frqfile-chr $FREQ \
    #         --print-coefficients \
    #         --print-delete-vals \
    #         --out $OUT
        
    #     Rscript $LDSC/ContinuousAnnotations/quantile_h2g.r \
    #         $BASELINE_QM \
    #         $OUT \
    #         $OUT.q5.txt 
    # fi
    
    # Run annot vs. ocr
    # echo "annot vs. ocr"
    # OUT="annotations/$NAME/out_ocr/$trait.$NAME"
    # if [ -f "$OUT.q5.txt" ]; then
    #     echo "Skipping because $OUT.q5.txt already exists..."
    # else
    #     python $LDSC/ldsc.py \
    #         --h2 $sumstat \
    #         --ref-ld-chr $LD_SCORE/$NAME.,$OCR_LDSCORE \
    #         --w-ld-chr $WEIGHTS \
    #         --overlap-annot \
    #         --frqfile-chr $FREQ \
    #         --print-coefficients \
    #         --print-delete-vals \
    #         --out $OUT 
        
    #     Rscript $LDSC/ContinuousAnnotations/quantile_h2g.r \
    #         $OCR_QM \
    #         $OUT \
    #         $OUT.q5.txt
    # fi

    # Run annot vs. ocr vs. baseline
    echo "annot vs. ocr vs. baseline"
    OUT="annotations/$NAME/out_ocr_baseline/$trait.$NAME"
    if [ -f "$OUT.q5.txt" ]; then
        echo "Skipping because $OUT.q5.txt already exists..."
    else
        python $LDSC/ldsc.py \
            --h2 $sumstat \
            --ref-ld-chr $LD_SCORE/$NAME.,$OCR_LDSCORE,$BASELINE \
            --w-ld-chr $WEIGHTS \
            --overlap-annot \
            --frqfile-chr $FREQ \
            --print-coefficients \
            --print-delete-vals \
            --out $OUT

        Rscript $LDSC/ContinuousAnnotations/quantile_h2g.r \
            $OCR_BASELINE_QM \
            $OUT \
            $OUT.q5.txt
    fi

    
done