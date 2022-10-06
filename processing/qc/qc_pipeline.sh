#!/usr/bin/env bash

while getopts i:c:o:e: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
        c) COHORT=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        e) ENVIRONMENT=${OPTARG};;
    esac
done

source /data/$USER/conda/etc/profile.d/conda.sh && conda activate $ENVIRONMENT

python3 /insert_wd_path/processing/qc/qc_pipeline.py --geno $INPUT --cohort $COHORT --out $OUTPUT
