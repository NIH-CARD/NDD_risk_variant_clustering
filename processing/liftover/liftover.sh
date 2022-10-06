#!/usr/bin/env bash

while getopts i:d:o: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
        d) OUTDIR=${OPTARG};;
        o) OUTPUT=${OPTARG};;
    esac
done

module load python/2.7
module load plink/1.9

plink --bfile $INPUT --recode --out $OUTDIR/recode

python liftOverPlink.py --map $OUTDIR/recode.map --out $OUTDIR/lift --chain insert_chain_file_path --bin insert_liftOver_binary_path

python rmBadLifts.py --map $OUTDIR/lift.map --out $OUTDIR/good_lifted.map --log $OUTDIR/bad_lifted.dat

module load python/3.7

awk '{print $2}' $OUTDIR/good_lifted.map > $OUTDIR/snplist.txt

plink --bfile $INPUT --extract $OUTDIR/snplist.txt --recode --out $OUTDIR/lifted

plink --ped $OUTDIR/lifted.ped --map $OUTDIR/good_lifted.map --make-bed --out $OUTPUT
