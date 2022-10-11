# NDD_risk_variant_clustering
Analysis of GWAS-significant SNPs across AD, PD, ALS, FTD and LBD

## Processing
### 1. Preprocessing
- ROSMAP/MayoRNAseq/MSBB: Concatenate individual chromosome vcf files and convert to PLINK binary files.
- ADNI: Merge individual chromosome PLINK binary files.
- DementiaSeq: Separate FTD, LBD, and ALS PLINK binary files.
- AMP-PD/ADSP: No processing needed to get PLINK binary files.
- All cohorts: Attach phenotype and sex data to .fam files.

### 2. QC
- Individual cohort QC run through GenoTools (https://github.com/dvitale199/GenoTools).
- ROSMAP/MayoRNAseq/MSBB, ADNI, AMP-PD are the cohorts that undergo genetic sex confirmation. Other cohorts did not have high-quality 23rd chromosome data available.

### 3. LiftOver
- From hg19 to hg38.
- ROSMAP/MayoRNAseq/MSBB, AD stats, PD stats, ALS stats and FTD stats.
- Using code adapted from https://genome.sph.umich.edu/wiki/LiftOver#Lift_genome_positions.

### 4. Merge
- Merge cohorts across common SNPs.
- Make phenotype/cohort information file.
- Use annovar to convert from chr:bp:ref:alt ro rsID (https://annovar.openbioinformatics.org/en/latest/).
- Run relatedness and duplicate prune as well as ancestry prune to ensure all samples are of EUR ancestry (again using GenoTools).
- Extract GWAS significance (i.e. p < 5e-08) variants from summary stats.
- Annotate gene-associated SNPs. When a SNP did not have a gene association, annotate the nearest gene.
- Sample 1000 cases from each disease.

### 5. Adjustment
- GenoML munging and PC adjustment (https://github.com/GenoML/genoml2).


## Analysis
### 1. Get UMAP Parameters
- Runs via swarm. Returns a .txt file listing the top parameter combinations based on the Disease ~ Cluster regressions.

### 2. Cluster Consistency
- Runs via swarm. Gets iterative cluster membership for multi-disease and disease-specific clusters

### 3. Analysis
- Queues the swarm jobs for Get UMAP Parameters and Cluster Consistency.
- Run full analysis. This needs to be separately run with disease set to None for multi-disease analysis and disease set to 'ad', 'pd', 'als', 'ftd', 'lbd' for disease-specific analyses.
- Note that for single disease analysis disease_regression() will fail.
