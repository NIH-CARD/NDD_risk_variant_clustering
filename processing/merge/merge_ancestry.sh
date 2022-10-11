#!/usr/bin/env bash

while getopts g:r:l:m:o:e: flag
do
    case "${flag}" in
        g) GENO=${OPTARG};;
        r) REF=${OPTARG};;
        l) LABELS=${OPTARG};;
        m) MODEL=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        e) ENVIRONMENT=${OPTARG};;
    esac
done

source /insert_conda_env_path/conda.sh && conda activate $ENVIRONMENT

python3 merge_ancestry.py --geno $GENO --ref $REF --ref_labels $LABELS --model $MODEL --out $OUTPUT
