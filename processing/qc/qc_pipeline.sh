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

source /insert_conda_env_path/conda.sh && conda activate $ENVIRONMENT

python3 /insert_wd_path/processing/qc/qc_pipeline.py --geno $INPUT --cohort $COHORT --out $OUTPUT
