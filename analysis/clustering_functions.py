import pandas as pd
import subprocess
import sys
import numpy as np
import os
import shutil
import joblib
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import MeanShift, DBSCAN
from sklearn.metrics import roc_auc_score, r2_score, mean_squared_error, accuracy_score
from sklearn.linear_model import LogisticRegression, LinearRegression 
from sklearn.ensemble import RandomForestRegressor
from umap import UMAP
import plotly.express as px
import plotly.io as pio
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import logit, ols, glm
import umap
import shap



def set_up(cases, disease=None):
    # isolate a single disease if requested
    if disease:
        cases = cases[cases['PHENO'] == disease]
        
    # get X data
    X = cases.drop(columns=['PHENO'], axis=1)
    
    # get y data
    y = cases['PHENO']
    
    return X, y



def umap_transform(X_train, X_test, y_cases, a, b, seed=None, cohorts=None):
    wd = 'insert_wd_path'
    
    # umap transform train and test data
    umap = UMAP(n_components=3, a=a, b=b, random_state=seed)
    X_train_umap = umap.fit_transform(X_train)
    X_test_umap = umap.transform(X_test)
    
    # get full umap data
    X_cases_umap = np.append(X_train_umap, X_test_umap, axis=0)
    
    return X_train_umap, X_test_umap



def fit_cluster(cluster, X_train, X_test, y_train, y_test, data_type, fname):
    # fit cluster and predict train and test clusters
    cluster.fit(X_train)
    clusters_train = cluster.predict(X_train)
    clusters_test = cluster.predict(X_test)
    print(pd.Series(clusters_train).value_counts())
    print(pd.Series(clusters_test).value_counts())
    
    # get full data
    cases = np.append(X_train[:,:3], X_test[:,:3], axis=0)
    phenotypes = np.append(y_train, y_test)
    clusters = np.append(clusters_train, clusters_test)
    
    # train dataframe
    train = pd.DataFrame(X_train)
    train['cluster'] = clusters_train
    train = train.reset_index(drop=True)
    y_train = y_train.reset_index(drop=True)
    train['pheno'] = y_train
    
    # test dataframe
    test = pd.DataFrame(X_test)
    test['cluster'] = clusters_test
    test = test.reset_index(drop=True)
    y_test = y_test.reset_index(drop=True)
    test['pheno'] = y_test
    
    return train, test



def expand_cluster_data(cluster_data):
    # removing cases that were not clustered (not an issue with MeanShift)
    cluster_data = cluster_data[cluster_data['cluster'] != -1]
    
    # get dummies for pheno and cluster membership for regression
    cluster_data = pd.concat([cluster_data, pd.get_dummies(cluster_data['pheno']), pd.get_dummies(cluster_data['cluster'], prefix='cluster')], axis=1)
    rename_dict = {}

    # get comp columns to rename
    for column in cluster_data.columns:
        if type(column) is int:
            rename_dict[column] = 'COMP' + str(column+1)

    # rename columns
    cluster_data = cluster_data.rename(columns=rename_dict)
    
    return cluster_data



def get_full_data(X_train, X_test, train_clusters, test_clusters):
    X_train = X_train.reset_index(drop=True)
    train_full = pd.concat([X_train, train_clusters], axis=1)
    
    X_test = X_test.reset_index(drop=True)
    test_full = pd.concat([X_test, test_clusters], axis=1)
    
    return train_full, test_full


    
def get_cluster_counts(data, fname):
    cols = ['Disease','Overall']
    cluster_cols = []
    for col in data.columns:
        if 'cluster_' in col:
            num = col.split('_')[1]
            cols.append(f'C{num}')
            cluster_cols.append(col)
    
    disease_cluster_representation = pd.DataFrame(columns=cols)
    disease_cluster_representation['Disease'] = ['ad','pd','lbd','ftd','als']

    data_dict = {}

    data_dict['Overall'] = dict(round(data['pheno'].value_counts(normalize=True), 5))
    
    for i in range(len(cluster_cols)):
        data_cluster = data[data[cluster_cols[i]] == 1]
        data_dict[f'C{i}'] = dict(round(data_cluster['pheno'].value_counts(normalize=True), 5))

    for cluster in data_dict:
        col_data = []
        for disease in disease_cluster_representation['Disease']:
            col_data.append(data_dict[cluster][disease])
        disease_cluster_representation[cluster] = col_data

    print(disease_cluster_representation)
    disease_cluster_representation.to_csv(f'results/{fname}', sep=',', index=False)
    
    

def check_cluster_membership(X_train, X_test, y_train, y_test, y_cases, y_ids, a, b, cluster, iterations, adjusted=None, prs=False):
    # set up data storage
    num_clusters = []
    cluster_membership_df = []
    cluster_membership_df = pd.DataFrame()
    cluster_membership_df['ID'] = y_ids
    
    # loop through passed num iterations
    for i in range(iterations):
        # umap transform
        X_train_umap, X_test_umap = umap_transform(X_train, X_test, y_cases, a, b, seed=None, cohorts=None)
        
        # fit clusters
        train_clusters, test_clusters = fit_cluster(cluster, X_train_umap, X_test_umap, y_train, y_test, 'UMAP', f'figures/umap_clusters_{i+1}.html')
    
        # get num clusters
        num_clusters.append(len(pd.Series(train_clusters['cluster']).value_counts()))
        
        # put cluster data into df
        cluster_data = np.append(train_clusters['cluster'], test_clusters['cluster'])
        cluster_membership_df[i+1] = cluster_data
    
    return cluster_membership_df, num_clusters



def disease_regression(train_data, test_data, data_type, standardize=False, simul=False):
    train_data_copy = train_data.copy(deep=True)
    test_data_copy = test_data.copy(deep=True)
    
    diseases = ['ad','pd','als','ftd','lbd']
    
    cluster_cols = []
    for column in train_data_copy.columns:
        if 'cluster_' in column:
            cluster_cols.append(column)
            
    comp_cols = []
    for column in train_data_copy.columns:
        if 'COMP' in column:
            comp_cols.append(column)
    
    if data_type == 'clusters':
        cols = cluster_cols
        axis = 'Cluster Membership'
    else:
        cols = comp_cols
        if standardize:
            scaler = StandardScaler()
            train_data_copy[comp_cols] = scaler.fit_transform(train_data_copy[comp_cols])
            test_data_copy[comp_cols] = scaler.fit_transform(test_data_copy[comp_cols])
            axis = 'Standardized Component'
        else:
            axis = 'Component'
    
    auc_data = []
    
    for disease in diseases:
        if simul:
            formula_str = ''
        
            for i in range(len(cols)):
                if data_type == 'clusters':
                    if cols[i] != cols[-2]:
                        if cols[i] != cols[-1]:
                            formula_str += str(cols[i]) + ' + '
                        else:
                            formula_str += str(cols[i])
                else:
                    if cols[i] != cols[-1]:
                            formula_str += str(cols[i]) + ' + '
                    else:
                        formula_str += str(cols[i])
            
            formula = (f'{disease} ~ {formula_str}')
            model = logit(formula=formula, data=train_data_copy).fit(disp=0)
            pred = model.predict(test_data_copy[cols])
            auc_data.append(roc_auc_score(test_data_copy[disease], pred))
        
        else:
            for col in cols:
                try:
                    formula = (f'{disease} ~ {col}')
                    model = logit(formula=formula, data=train_data_copy).fit(disp=0)
                    pred = model.predict(test_data_copy[col])
                    auc_data.append(roc_auc_score(test_data_copy[disease], pred))  
                except (np.linalg.LinAlgError, KeyError) as e:
                    auc_data.append(0)
    
    print('Average AUC')
    print(np.mean(auc_data))
    
    return auc_data
