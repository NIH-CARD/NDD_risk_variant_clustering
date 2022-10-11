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

from clustering_functions import set_up, fit_cluster, expand_cluster_data, disease_regression

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
    pheno.columns = ['FID','ID','PHENO','COHORT']

    cases_controls = adjusted.merge(pheno[['ID','PHENO']], how='inner', on=['ID'])
    cases = cases_controls[cases_controls['PHENO'] != 'control'].reset_index(drop=True)

    X, y = set_up(cases)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=None)

    X_train = X_train.drop(columns=['ID'], axis=1)
    X_test = X_test.drop(columns=['ID'], axis=1)

    y_cases = np.append(y_train, y_test)

    a_vals = np.arange(0.25, 3.25, 0.25)
    b_vals = np.arange(0.25, 3.25, 0.25)

    cluster_auc_dict = {}

    for a in a_vals:
        for b in b_vals:

            print(f'{a}, {b}')
            umap = UMAP(n_components=3, a=a, b=b)
            X_train_umap = umap.fit_transform(X_train)
            X_test_umap = umap.transform(X_test)

            cluster = MeanShift()
            train_data, test_data = fit_cluster(cluster, X_train_umap, X_test_umap, y_train, y_test, None, None)

            train_data = expand_cluster_data(train_data)
            test_data = expand_cluster_data(test_data)

            cluster_auc = disease_regression(train_data, test_data, 'clusters')

            key = f'a: {str(a)}, b: {str(b)}'

            cluster_auc_dict[key] = np.mean(cluster_auc)

    for key in cluster_auc_dict:
        if cluster_auc_dict[key] == np.nan:
            cluster_auc_dict = cluster_auc_dict.pop(key)

    print('CLUSTERS')
    print(cluster_auc_dict)
    print()
    sorted_clusters = sorted(cluster_auc_dict, key=cluster_auc_dict.get, reverse=True)[:10]
    print(sorted_clusters)
    print()

    with open(f'params/umap_params_downsampled_adj_random_state_None.txt', 'w') as f:
        f.write('\nCLUSTERS\n')
        for key in sorted_clusters:
            f.write(f'{key}: {cluster_auc_dict[key]}\n')
        f.close()



