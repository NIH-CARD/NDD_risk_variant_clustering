#!/usr/bin/env bash

while getopts i:c:p:o: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
        c) CHROM=${OPTARG};;
        p) POS=${OPTARG};;
        o) OUT=${OPTARG};;
    esac
done

export R_LIBS_USER="/data/$USER/R/4.1"

module load R/4.1

Rscript --vanilla /insert_wd_path/processing/liftover/liftover_summary_stats.R --summary_stats $INPUT --from_build "hg19" --to_build "hg38" --chr "$CHROM" --pos "$POS" --out $OUT