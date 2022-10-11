import pandas as pd
import argparse
import shutil
import os

# local imports
from QC.qc import callrate_prune, het_prune, sex_prune, related_prune, variant_prune, miss_rates
from QC.utils import shell_do


parser = argparse.ArgumentParser(description='Arguments for Genotyping QC (data in Plink .bim/.bam/.fam format)')
parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
parser.add_argument('--cohort', type=str, default='nope', help='Cohort being QC\'d')
parser.add_argument('--out', type=str, default='nope', help='Prefix for output (including path)')

args = parser.parse_args()

geno_path = args.geno
cohort = args.cohort
out_path = args.out

# making qc dir if it doesn't already exist
dir_name = os.path.dirname(geno_path)
qc_dir = f'{dir_name}/qc'
os.makedirs(qc_dir, exist_ok=True)

# output list
outputs = []

# metrics path
metrics_path = f'{out_path}_metrics.txt'

missing_path = f'{out_path}_missing'
avg_miss = miss_rates(geno_path, missing_path)

# writing missing rates to metrics file
with open(metrics_path, 'a') as f:
    for key, value in avg_miss.items():
        f.write(f'{key}: {value}\n')
    f.close()

callrate_out = f'{out_path}_callrate'
callrate = callrate_prune(geno_path, callrate_out, mind=0.05)
outputs.append(callrate)

if cohort not in ['FTD','LBD','ALS','ADSP']:
    sex_out = f'{callrate_out}_sex'
    sex = sex_prune(callrate_out, sex_out)
    outputs.append(sex)
    related_out = f'{sex_out}_related'
    related = related_prune(sex_out, related_out)
else:
    related_out = f'{callrate_out}_related'
    related = related_prune(callrate_out, related_out)

outputs.append(related)

het_out = f'{related_out}_het'
het = het_prune(related_out, het_out)
outputs.append(het)

if het['pass']:
    variant = variant_prune(het_out, out_path)
else:
    variant = variant_prune(related_out, out_path)
outputs.append(variant)

# writing metrics from other steps to metrics file
with open(metrics_path, 'a') as f:
    for out in outputs:
        f.write(f"{out['step']}\n")
        for key, value in out['metrics'].items():
            f.write(f'\t{key}: {value}\n')
    f.close()