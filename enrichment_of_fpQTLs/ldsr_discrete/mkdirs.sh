#!/bin/bash

for ANNOT in TOBIAS TOBIAS_cpm_normalized; do

    echo "Making dirs for $ANNOT..."
    mkdir annotations/$ANNOT
    mkdir annotations/$ANNOT/out
    mkdir annotations/$ANNOT/ld_score

    echo "Making dirs for ${ANNOT}_covariates..."
    mkdir annotations/${ANNOT}_covariates
    mkdir annotations/${ANNOT}_covariates/out
    mkdir annotations/${ANNOT}_covariates/ld_score

done