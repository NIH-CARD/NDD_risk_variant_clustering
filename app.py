#upload some libraries
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import seaborn as sns
from PIL import Image
import datetime

def plot_3d(labeled_df, color, symbol=None, x='UMAP1', y='UMAP2', z='UMAP3', title=None, x_range=None, y_range=None, z_range=None):
    '''
    Parameters: 
    labeled_df (Pandas dataframe): labeled ancestry dataframe
    color (string): color of ancestry label. column name containing labels for ancestry in labeled_pcs_df
    symbol (string): symbol of secondary label (for example, predicted vs reference ancestry). default: None
    plot_out (string): filename to output filename for .png and .html plotly images
    x (string): column name of x-dimension
    y (string): column name of y-dimension
    z (string): column name of z-dimension
    title (string, optional): title of output scatterplot
    x_range (list of floats [min, max], optional): range for x-axis
    y_range (list of floats [min, max], optional): range for y-axis
    z_range (list of floats [min, max], optional): range for z-axis

    Returns:
    3-D scatterplot (plotly.express.scatter_3d). If plot_out included, will write .png static image and .html interactive to plot_out filename
        
    '''    
    fig = px.scatter_3d(
        labeled_df,
        x=x,
        y=y,
        z=z,
        color=color,
        symbol=symbol,
        title=title,
        color_discrete_sequence=px.colors.qualitative.Bold,
        range_x=x_range,
        range_y=y_range,
        range_z=z_range
    )

    fig.update_traces(marker={'size': 3})

    st.plotly_chart(fig)

def plot_prs_distributions_by_disease(full_data):
    # isolate each cluster
    cluster_0 = full_data[full_data['cluster'] == 0]
    cluster_1 = full_data[full_data['cluster'] == 1]
    cluster_2 = full_data[full_data['cluster'] == 2]
    
    diseases = ['ad','pd','als','lbd','ftd']
    
    fig, axs = plt.subplots(5)
    
    plt.subplots_adjust(hspace=0.35)
    
    fig.set_figwidth(6)
    fig.set_figheight(18)
    
    i = 0
    
    for disease in diseases:
        sns.kdeplot(data=full_data, x=f'PRS_STD_{disease}', hue='cluster', palette='bright', common_norm=False, ax=axs[i])
        axs[i].set_xlabel(f'{disease.upper()} PRS')
        axs[i].set_title(f'{disease.upper()} PRS Distributions')
        axs[i].axvline(x=0, c='black')
        i += 1
    
    
    st.pyplot(fig)

st.set_page_config(
    layout = 'wide'
)

st.markdown("""
    <style>
    .big-font {
        font-family:Helvetica; color:#0f557a; font-size:55px !important;
    }
    </style>
    """, unsafe_allow_html=True)

st.markdown("""
    <style>
    .small-font {
        font-family:Helvetica; color:#0f557a; font-size:20px !important;
    }
    </style>
    """, unsafe_allow_html=True)

st.markdown("""
    <style>
    .tiny-font {
        font-family:Helvetica; color:#0f557a; font-size:12px !important;
    }
    </style>
    """, unsafe_allow_html=True)

##Background color
def local_css(file_name):
    with open(file_name) as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
local_css(f'streamlit/style/style.css')
####################### HEAD ##############################################

title, head_1 = st.columns([1, 0.2])

with title:
    st.markdown('<p class="big-font">NDD Risk Variant Clustering</p>', unsafe_allow_html=True)

card = Image.open(f'streamlit/style/card.jpeg')
head_1.image(card, width=120)

########################  SIDE BAR #########################################
st.sidebar.markdown('**Stats/Graph Selection**', unsafe_allow_html=True)
selected_metrics = st.sidebar.selectbox(label="Disease Selection", options=['Multi-disease','AD','PD','ALS','FTD','LBD'])

if selected_metrics == 'Multi-disease':
    data_path = 'streamlit/data/full.csv'

else:
    data_path = f'streamlit/data/{selected_metrics.lower()}.csv'

data = pd.read_csv(data_path, sep=',')
data['cluster'] = data['cluster'].astype(str)

if selected_metrics == 'Multi-disease':
    st.markdown(f'<p class="small-font">{selected_metrics} UMAP Embedding Labeled by Disease</p>', unsafe_allow_html=True)
    plot_3d(data, 'pheno')

    st.markdown(f'<p class="small-font">{selected_metrics} Cluster Membership</p>', unsafe_allow_html=True)
    plot_3d(data, 'cluster')

else:
    st.markdown(f'<p class="small-font">{selected_metrics} Cluster Membership</p>', unsafe_allow_html=True)
    plot_3d(data,'cluster')

tables, prs_plot = st.columns([1, 0.5])

with prs_plot:
    st.markdown('<p class="small-font">PRS Distributions by Cluster</p>', unsafe_allow_html=True)
    plot_prs_distributions_by_disease(data)

with tables:

    if selected_metrics == 'Multi-disease':
        st.markdown('<p class="small-font">Disease association summary stats per cluster/Frequency per cluster</p>', unsafe_allow_html=True)
        cluster_regression_path = 'streamlit/data/full_disease_regression.csv'
        cluster_regression = pd.read_csv(cluster_regression_path)
        cluster_regression['OR'] = cluster_regression['OR'].map('{:.3f}'.format)
        cluster_regression['BETA'] = cluster_regression['BETA'].map('{:.3f}'.format)
        cluster_regression['SE'] = cluster_regression['SE'].map('{:.3f}'.format)
        st.dataframe(cluster_regression, height=384)
        st.markdown('<p class="tiny-font">* denotes a p-value < 0.05 for the frequency increase/decrease in a certain disease \
            status per cluster compared to the null estimate of 20%</p>', unsafe_allow_html=True)

        st.markdown('<p class="small-font">PRS association summary stats per cluster</p>', unsafe_allow_html=True)
        prs_regression_path = 'streamlit/data/full_prs_regression.csv'
        prs_regression = pd.read_csv(prs_regression_path)
        prs_regression['OR'] = prs_regression['OR'].map('{:.3f}'.format)
        prs_regression['BETA'] = prs_regression['BETA'].map('{:.3f}'.format)
        prs_regression['SE'] = prs_regression['SE'].map('{:.3f}'.format)
        st.dataframe(prs_regression, height=384)

    else:
        st.markdown(f'<p class="small-font">{selected_metrics} cluster counts</p>', unsafe_allow_html=True)
        counts_path = f'streamlit/data/{selected_metrics.lower()}_cluster_counts' 
        counts = pd.read_csv(f'{counts_path}.csv')
        st.dataframe(counts)

        st.markdown(f'<p class="small-font">{selected_metrics} PRS regression summary stats per cluster</p>', unsafe_allow_html=True)
        prs_regression_path = f'streamlit/data/{selected_metrics.lower()}_prs_regression' 
        prs_regression = pd.read_csv(f'{prs_regression_path}.csv')
        prs_regression['OR'] = prs_regression['OR'].map('{:.3f}'.format)
        prs_regression['BETA'] = prs_regression['BETA'].map('{:.3f}'.format)
        prs_regression['SE'] = prs_regression['SE'].map('{:.3f}'.format)
        st.dataframe(prs_regression)

        st.markdown(f'<p class="small-font">Cross-disease PRS associations [mean (sd)]</p>', unsafe_allow_html=True)
        prs_assoc_path = f'streamlit/data/{selected_metrics.lower()}_prs_associations' 
        prs_assoc = pd.read_csv(f'{prs_assoc_path}.csv')
        st.dataframe(prs_assoc)
        st.markdown('<p class="tiny-font">* denotes a p-value < 0.05 for the deviation from the distribution (mean=0, sd=1) \
            that the PRS were standardized to', unsafe_allow_html=True)

