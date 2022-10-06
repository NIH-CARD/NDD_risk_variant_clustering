import pandas as pd
import argparse
import shutil
import os

from Ancestry.ancestry import run_ancestry

parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--ref', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format reference genotype file, everything before the *.bed/bim/fam.')
parser.add_argument('--ref_labels', type=str, default='nope', help='tab-separated plink-style IDs with ancestry label (FID  IID label) with no header')
parser.add_argument('--model', type=str, default=None, help='Path to pickle file with trained ancestry model for passed reference panel')
parser.add_argument('--out', type=str, default='nope', help='Prefix for output (including path)')

args= parser.parse_args()

geno_path = args.geno
ref_panel = args.ref
ref_labels = args.ref_labels
model_path = args.model
out_path = args.out

# making ancestry dir if it doesn't already exist
dir_name = os.path.dirname(geno_path)
ancestry_dir = f'{dir_name}/ancestry'
os.makedirs(ancestry_dir, exist_ok=True)

if os.path.isfile(model_path):
    ancestry = run_ancestry(geno_path=geno_path, out_path=out_path, ref_panel=ref_panel, ref_labels=ref_labels, model_path=model_path)
else:
    ancestry = run_ancestry(geno_path=geno_path, out_path=out_path, ref_panel=ref_panel, ref_labels=ref_labels, model_path=None)