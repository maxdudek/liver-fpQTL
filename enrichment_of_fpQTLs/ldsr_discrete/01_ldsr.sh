#!/bin/bash
#SBATCH --job-name="ldsr_fpQTLs"
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 3-00:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out



source ~/.bashrc
for i in $(seq ${CONDA_SHLVL}); do
    conda deactivate
done
conda activate /mnt/isilon/sfgi/programs/miniconda3/envs/ldsc

module load bzip2

LDSC="/mnt/isilon/sfgi/programs/ldsc"

# Test LDSC environment
which python
#python $LDSC/make_annot.py -h
#python $LDSC/ldsc.py -h
#exit

# Phase 3 - hg38
# PHASE3="/mnt/isilon/sfgi/pahlm/annotationFiles/LDSR_1k_phase3_ref"
# BFILE="$PHASE3/hg38_liftover/plink_files/1000G.EUR.hg38"
# PRINT_SNPS="$PHASE3/hapmap3_snps/hm"
# BASELINE="$PHASE3/hg38_liftover/baselineLD_v2.2/baselineLD."
# WEIGHTS="$PHASE3/hg38_liftover/weights/weights.hm3_noMHC."
# FREQ="$PHASE3/hg38_liftover/plink_files/1000G.EUR.hg38."

# Khanh hg38
cdir=/mnt/isilon/sfgi/trangk/analyses/grant/ldsc
BFILE=$cdir/plink_files_hg38/1000G.EUR.hg38
FREQ=$cdir/plink_files_hg38/1000G.EUR.hg38.
BASELINE=$cdir/baselineLD_v2.2_hg38/baselineLD.
WEIGHTS=$cdir/weights_hg38/weights.hm3_noMHC.
PRINT_SNPS=$cdir/list.txt

# echo "plink_files_hg38/1000G.EUR.hg38.1.bim:"
# wc -l $BFILE.1.bim
# echo ""

# echo "baselineLD_v2.2_hg38/baselineLD.1.annot.gz:"
# gunzip < "${BASELINE}1.annot.gz" | wc -l 
# echo ""

# echo "weights_hg38/weights.hm3_noMHC.1.l2.ldscore.gz:"
# gunzip < "${WEIGHTS}1.l2.ldscore.gz" | wc -l 
# echo ""

# echo "plink_files_hg38/1000G.EUR.hg38.1.frq:"
# wc -l "${FREQ}1.frq"
# echo ""

# Run to reset results
# for n in annotations/* ; do
#     arrIN=(${n//\// })
#     NAME=${arrIN[-1]}
#     
#     rm annotations/$NAME/out/*
# done

# for n in annotations/* ; do
#     arrIN=(${n//\// })
#     NAME=${arrIN[-1]}

    
# for NAME in HepG2_TOBIAS_footprints iPSC-hep_TOBIAS_footprints iPSC_TOBIAS_footprints brandon_liver_TOBIAS_footprints HepG2_ChIP_TFs hepatocyte_ChIP_TFs liver_ChIP_TFs brandon_liver_ocr ; do 
# for NAME in PRINT_no_gaussian PRINT_no_gaussian_FDR10; do
# for NAME in brandon_liver_ocr PRINT_beta PRINT_beta_FDR10 PRINT_no_gaussian PRINT_no_gaussian_FDR10; do
for NAME in PRINT_beta_cloglog PRINT_beta_cloglog_FDR10 PRINT_beta_cloglog_phi PRINT_beta_cloglog_phi_FDR1 PRINT_beta_logit PRINT_beta_logit_FDR10 PRINT_beta_probit PRINT_beta_probit_FDR10; do
# for n in annotations/PRINT_beta_* ; do
#     arrIN=(${n//\// })
#     NAME=${arrIN[-1]}

    BED_IN="annotations/$NAME/$NAME.bed"
    LD_SCORE="annotations/$NAME/ld_score"

    mkdir -p annotations/$NAME/ld_score
    mkdir -p annotations/$NAME/out_baseline
    mkdir -p annotations/$NAME/out_ocr
    mkdir -p annotations/$NAME/out_ocr_baseline

    echo "Running LDSR for functional annotation '$NAME'"

    # Create annotation files
    if [ -f "$LD_SCORE/$NAME.22.l2.ldscore.gz" ]; then
        echo "Skipping annotation creation for $NAME because $LD_SCORE/$NAME.22.l2.ldscore.gz already exists..." 
    else
        # Bed to Annot
        for j in {1..22}; do
        python $LDSC/make_annot.py \
            --bed-file $BED_IN \
            --bimfile $BFILE.$j.bim \
            --annot-file $LD_SCORE/$NAME.$j.annot.gz
        done

        #Annot to LD files
        for j in {1..22}; do
        python $LDSC/ldsc.py \
            --l2 \
            --bfile $BFILE.$j \
            --ld-wind-cm 1 \
            --annot $LD_SCORE/$NAME.$j.annot.gz \
            --thin-annot \
            --out $LD_SCORE/$NAME.$j \
            --print-snps $PRINT_SNPS
        done
    fi

    # echo "annotations/$NAME/ld_score/TOBIAS.1.annot.gz:"
    # gunzip < "$LD_SCORE/$NAME.1.annot.gz" | wc -l
    # echo ""

    # echo "annotations/$NAME/ld_score/TOBIAS.1.l2.ldscore.gz:"
    # gunzip < "$LD_SCORE/$NAME.1.l2.ldscore.gz" | wc -l
    # echo ""


    # Run LDSR, liver curated sumstats
    SUMSTATS="/mnt/isilon/sfgi/dudekm/raw_data/GWAS/sumstats/munged_sumstats"
    for sumstat in "$SUMSTATS"/*.sumstats.gz; do
        echo "Running regression on $sumstat..."

        # */NAFLD.sumstats.gz --> NAFLD
        arrIN=(${sumstat//\// })
        arrIN=(${arrIN[-1]//.sumstats./ })
        trait=${arrIN[0]}

        # Run annot vs. baseline
        echo "annot vs. baseline"
        OUT="annotations/$NAME/out_baseline/$trait.$NAME"
        if [ -f "$OUT.results" ]; then
            echo "Skipping because $OUT.results already exists..." 
        else
            python $LDSC/ldsc.py \
                --h2 $sumstat \
                --ref-ld-chr $LD_SCORE/$NAME.,$BASELINE \
                --w-ld-chr $WEIGHTS \
                --overlap-annot \
                --frqfile-chr $FREQ \
                --print-coefficients \
                --out $OUT
        fi

        if [[ "$NAME" == "brandon_liver_ocr" ]]; then 
            continue
        fi
        
        # Run annot vs. ocr
        OCR_LDSCORE="annotations/brandon_liver_ocr/ld_score/brandon_liver_ocr."
        echo "annot vs. ocr"
        OUT="annotations/$NAME/out_ocr/$trait.$NAME"
        if [ -f "$OUT.results" ]; then
            echo "Skipping because $OUT.results already exists..." 
        else
            python $LDSC/ldsc.py \
                --h2 $sumstat \
                --ref-ld-chr $LD_SCORE/$NAME.,$OCR_LDSCORE \
                --w-ld-chr $WEIGHTS \
                --overlap-annot \
                --frqfile-chr $FREQ \
                --print-coefficients \
                --out $OUT
        fi

        # Run annot vs. ocr vs. baseline
        echo "annot vs. ocr vs. baseline"
        OUT="annotations/$NAME/out_ocr_baseline/$trait.$NAME"
        if [ -f "$OUT.results" ]; then
            echo "Skipping because $OUT.results already exists..." 
        else
            python $LDSC/ldsc.py \
                --h2 $sumstat \
                --ref-ld-chr $LD_SCORE/$NAME.,$OCR_LDSCORE,$BASELINE \
                --w-ld-chr $WEIGHTS \
                --overlap-annot \
                --frqfile-chr $FREQ \
                --print-coefficients \
                --out $OUT
        fi
    done

    # Alkes Price group sumstats
    # SUMSTATS="/mnt/isilon/sfgi/dudekm/raw_data/GWAS/sumstats/alkesgroup_sumstats"
    # for sumstat in "$SUMSTATS"/*.sumstats.gz; do
    #     echo "Running regression on $sumstat..."

    #     # */NAFLD.sumstats --> NAFLD
    #     arrIN=(${sumstat//\// })
    #     arrIN=(${arrIN[-1]//.sumstats/ })
    #     trait=${arrIN[0]}

    #     OUT="annotations/$NAME/out/$trait.$NAME"

    #     if [ -f "$OUT.results" ]; then
    #         echo "Skipping because $OUT.results already exists..."
    #         continue 
    #     fi

    #     # #run ldscore regression
    #     # #h2 outpt from munge sumstats
    #     python $LDSC/ldsc.py \
    #         --h2 $sumstat \
    #         --ref-ld-chr $LD_SCORE/$NAME.,$BASELINE \
    #         --w-ld-chr $WEIGHTS \
    #         --overlap-annot \
    #         --frqfile-chr $FREQ \
    #         --print-coefficients \
    #         --out $OUT
    # done
done

