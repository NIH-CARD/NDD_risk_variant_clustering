#!/usr/bin/env bash

while getopts i:o:e: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        e) ENVIRONMENT=${OPTARG};;
    esac
done

source /insert_conda_env_path/conda.sh && conda activate $ENVIRONMENT

python3 /insert_wd_path/processing/merge/merge_qc.py --geno $INPUT --out $OUTPUT
