# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 14:13:39 2022

@author: Youri LRM
"""

#%%

import os
print(os.getcwd())

#%%

import numpy as np
import pandas as pd
import scanpy as sc

#%%

# tuto
# https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
# TENSORFLOW
# HDF5

# Seurat:
# LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
# This is then natural-log transformed using log1p.
# x = log(((500 * 10000) / 50000 ) + 1)


wd = "data/Glimmunology_GBM_1/Glioma_Y_and_O/Levi2_Glioma_Y/outs/raw_feature_bc_matrix"


adata = sc.read_10x_mtx(
    wd,  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   

#%%

adata.var_names_make_unique()


sc.pl.highest_expr_genes(adata, n_top=20, )


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# The small it's set, the more it excludes?
#adata = adata[adata.obs.n_genes_by_counts < 1500, :]
#adata = adata[adata.obs.pct_counts_mt < 4, :]
adata = adata[adata.obs.n_genes_by_counts > 1500, ]
adata = adata[adata.obs.pct_counts_mt < 4,]


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')



sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


# Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells.
sc.pp.normalize_total(adata, target_sum=1e4)


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)


sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=5)






