import pandas as pd
import argparse
import shutil
import os

# local imports
from QC.qc import related_prune
from QC.utils import shell_do


parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--out', type=str, default='nope', help='Prefix for output (including path)')

args = parser.parse_args()

geno_path = args.geno
out_path = args.out

# making qc dir if it doesn't already exist
dir_name = os.path.dirname(geno_path)
qc_dir = f'{dir_name}/qc'
os.makedirs(qc_dir, exist_ok=True)

related = related_prune(geno_path, out_path)
print(related)