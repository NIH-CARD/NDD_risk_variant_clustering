import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.cluster import MeanShift
from sklearn.metrics import roc_auc_score
from umap import UMAP
import statsmodels.api as sm
from statsmodels.formula.api import logit
import statsmodels
import seaborn as sns
import sklearn
import numba
import umap
import sys

from clustering_functions import set_up, check_cluster_membership

if __name__ == '__main__':

    print(umap.__version__)
    print(np.__version__)
    print(sklearn.__version__)
    print(numba.__version__)
    print(sys.executable)
    print(sys.version)
    print(sys.version_info)

    wd = 'insert_path'
    adjusted_path = f'{wd}/processing/adjustment/downsampled/downsampled_gwas5e08_ADJUSTED10PCs.csv'
    pheno_path = f'{wd}/merged_genotypes/all_cohorts_phenotype.txt'

    adjusted = pd.read_csv(adjusted_path, sep=',')

    pheno = pd.read_csv(pheno_path, sep='\s+', header=None)
    pheno.columns = ['ID','IID','PHENO','COHORT']
    pheno = pheno.drop_duplicates(subset=['ID','IID'], ignore_index=True)

    cases_controls = adjusted.merge(pheno[['ID','PHENO','COHORT']], how='inner', on=['ID'])
    cases_controls = cases_controls.drop_duplicates(ignore_index=True)

    cases = cases_controls[cases_controls['PHENO'] != 'control'].reset_index(drop=True)

    diseases = [None,'ad','pd','als','ftd','lbd']
    
    for disease in diseases:

        X, y = set_up(cases, disease)
        col_names = X.columns

        for i in range(15):
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=None)
            train_id = X_train['ID']
            pd.concat([train_id, train_id], axis=1).to_csv('train_id.txt', sep='\t', index=False, header=False)
            test_id = X_test['ID']
            pd.concat([test_id, test_id], axis=1).to_csv('test_id.txt', sep='\t', index=False, header=False)
            train_cohort = X_train['COHORT']
            test_cohort = X_test['COHORT']
            X_train = X_train.drop(columns=['ID','COHORT'], axis=1)
            X_test = X_test.drop(columns=['ID','COHORT'], axis=1)
            y_cases = np.append(y_train, y_test)
            y_ids = np.append(train_id, test_id)
            y_cohorts = np.append(train_cohort, test_cohort)

            cluster = MeanShift()
            iterations = 100
            cluster_membership_df, num_clusters = check_cluster_membership(X_train, X_test, y_train, y_test, y_cases, y_ids, 2.75, 0.75, cluster, iterations)

            cluster_membership_df['pheno'] = y_cases
            print(cluster_membership_df.head())
            print(pd.Series(num_clusters).value_counts())

            if disease:
                cluster_membership_df.to_csv(f'bootstrap/cluster_membership_df_{str(iterations)}_rs_None_{str(i+1)}_{disease}.txt', sep='\t', index=False)
            else:
                cluster_membership_df.to_csv(f'bootstrap/cluster_membership_df_{str(iterations)}_rs_None_{str(i+1)}.txt', sep='\t', index=False)