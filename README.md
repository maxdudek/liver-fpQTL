# Liver fpQTL Discovery

This repository contains scripts used in the manuscript: \
[Characterization of non-coding variants associated with transcription factor binding through ATAC-seq-defined footprint QTLs in liver](https://doi.org/10.1101/2024.09.24.614730) (Preprint) 

Raw data files will be freely availible online pending publication.

Currently, all raw and intermediate data files are availible on reasonable request. Any questions can be sent to maxdudek@gmail.com.

## Directories

Below is the list of directories corresponding to the order of analyses presented in the paper. Most scripts within directories are numbered according to execution order. Code that corresponds to figures in the manuscript is labeled with a comment e.g. `#MANUSCRIPT - Figure 2A`

### 1. `PRINT`

Calculate footprint (FP) scores in all variants across all samples using ATAC-seq data.

### 2. `regression`

Normalized data and run regressions (and simulated regressions) comparing genotype to FP score to find fpQTLs.

### 3. `enrichment_of_fpQTLs`

Calculate the enrichment of fpQTLs in various functional annotations.

#### 3a. `enrichment_of_fpQTLs/functional_annotations`

Calculate the enrichment of fpQTLs in ENCODE annotations and TSS proximal regions

#### 3b. `enrichment_of_fpQTLs/chip-seq_overlap`

Label variants based on overlap with ChIP-seq peaks.

#### 3c. `enrichment_of_fpQTLs/ADASTRA`

Intersect variants with allele-specific binding data from ADASTRA.

#### 3d. `enrichment_of_fpQTLs/QTLs`

Label fpQTLs with eQTL and caQTL effect sizes, and overlap with caQTL credible sets.

#### 3e. `enrichment_of_fpQTLs/tables_2x2`

Label fpQTLs with GWAS significance and eQTL/caQTL significance, and construct 2x2 tables to calculate enrichment with GWAS proxies.

#### 3f. `enrichment_of_fpQTLs/ldsr_discrete`

Calculate the enrichment of trait heritability in fpQTLs using stratified linkage disequilibirum score regression (S-LDSR).

#### 3g. `enrichment_of_fpQTLs/ldsr_continuous`

Run S-LDSR using a continuous annotation based on *P*-values from fpQTL regresison

### 4. `allele_specific_binding`

Calculate allele-specific binding FP scores within heterozygous samples.

### 5. `motif_analysis`

Label variants based on overlap with known TF binding motifs and determine their concordance with fpQTL effect size.

### 6. `annotation`

Combine multiple fpQTL annotations in order to create the Table of significnat fpQTLs and to facilitate downstream analysis. Also generate gene track plots of GWAS loci as fine-mapping examples.

### 7. `cutsite_plotting`

Plot insertion pileups, heatmaps, and boxplots for all fpQTLs

### 8. `misc`

#### 8a. `misc/manhattan_peaks`

Generate the manhattan plot and zoom in to asses specific features.

#### 8b. `misc/chr17_inversion`

Analyze the effect of the chr17 inversion on fpQTL calling. Run an additional regression which accounts for the inversion haplotypes.

#### 8c. `misc/telemeric_fpQTLs`

Analyze fpQTLs in "telomeric-neighboring" regions (within 2 Mb of a chromosome end)
