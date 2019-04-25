#!/usr/bin/env python3

"""
> cluster_analysis.py <

Plots a correlation heatmap using seaborn.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.spatial.distance
import seaborn as sns

# define strains; underscore denotes union
strains = ['CC7', 'H2', 'RS', 'SSB01', 'CC7_SSB01', 'H2_SSB01', 'all']

# read in 16S abundances
data = pd.read_table('abundance.family.tsv', usecols=range(1,49))

for s in strains:
    # modify regex to catch multiple strains if needed. only use regex when
    # not 'all'
    if s != 'all':
        r = s.replace('_', '_|') + '_'
        data_subset = data.filter(regex=r)
    else:
        data_subset = data
    
    # remove water samples
    #data_subset = data_subset.filter(regex='^(?!H2O)')
    
    # drop rows with all zeroes
    data_subset = data_subset.loc[(data_subset!=0).any(axis=1)]
    
    # distance matrix: kendall or pearson
    #data_corr = data_subset.corr(method='kendall') 
    #data_corr = data_subset.corr(method='pearson')
    
    # distance matrix: braycurtis
    data_corr = scipy.spatial.distance.squareform(
        scipy.spatial.distance.pdist(data_subset.transpose(), 'braycurtis'))
    data_corr = pd.DataFrame(data_corr, index=data_subset.columns, columns=data_subset.columns)
    
    mask = np.zeros((len(data_corr), len(data_corr)), 'int8')
    np.fill_diagonal(mask, 1)
    
    # plot the heatmap!
    fig, ax = plt.subplots(figsize=(20, 20))
    
    # reduce font_scale for 'all'
    if s != 'all':
        sns.set(font_scale=1)
    else:
        sns.set(font_scale=0.8)
    
    cg = sns.clustermap(data_corr, method='ward')#,
                        vmin=0.35, vmax=0.85, mask=mask,
                        cmap='YlGnBu', annot=False, linewidth=0.5)
                        #cbar_kws={'ticks': [0.35, 0.45, 0.55, 0.65, 0.75]})
    
    # save figure
    fig.tight_layout()
    fig = plt.gcf()
    
    # without bbox_inches, the saved figure has truncated axes.
    output_filename = 'corr_plot.{}.braycurtis.pdf'.format(s)
    fig.savefig(output_filename, bbox_inches='tight')
    
    temp = data_corr - mask
    print (data_corr.min().min(), temp.max().max())
