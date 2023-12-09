import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from copy import copy

plt.rcParams["pdf.use14corefonts"] = True
sc.set_figure_params(vector_friendly=True)

##check the expression of common cytokines across decidual cell types using scRNA-seq data from Vento-Tormo et al., Nature, 2018
adata = sc.read('nature_2018_10X.h5ad')
adata = adata.raw.to_adata()
adata.obs['celltype'] = adata.obs['annotation'].astype(str)
ind = ((adata.obs['location'] != 'Blood') & (~adata.obs['celltype'].isin(['Granulocytes', 'MO', 'NK CD16+', 'NK CD16-', 'Plasma'])))
adata = adata[ind, :].copy()

ptys = ['dNK1', 'dNK2', 'dNK3', 'dNK p', 'dM1', 'dM2', 'dM3', 'DC1', 'DC2', 'dS1', 'dS2', 'dS3', 'dP1', 'dP2',
        'Endo (m)', 'Epi1', 'Epi2', 'ILC3', 'Tcells']
ind = adata.obs['celltype'].isin(ptys)
adata = adata[ind, :].copy()
sc.pp.filter_genes(adata, min_cells=20)

ggs = pd.read_csv('all_cytokines.txt', header=None)
ggs.index = ggs[0].values
ggs = ggs.index.intersection(adata.var_names)
print(len(ggs))

fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(25, 6), gridspec_kw=dict(bottom=0.1, top=0.99))
sc.pl.dotplot(adata, var_names=ggs, groupby='celltype', standard_scale='var', categories_order=ptys, cmap='GnBu', ax=ax, show=False)
fig.savefig('common_cytokines_decidual_celltypes.pdf')

##plot the 4 uNK-restricted cytokines
#This h5ad file was downloaded from https://data.humancellatlas.org/ to obtain the original UMAP coordinates
adata = sc.read('03ef13ff-c0ce-5f9f-a5c6-58cdeb288e86/decidua.h5ad')
adata.var_names = adata.var_names.str.replace('_ENSG.*', '')
adata.var_names_make_unique()
ind = ((adata.obs['Location'] != 'Blood') & (~np.isin(adata.obs['CellType'], ['Granulocytes', 'MO', 'NK CD16+', 'NK CD16-', 'Plasma'])))
adata = adata[ind, :].copy()
adata.obs['celltype'] = adata.obs['CellType'].astype(str)
adata.obs.loc[adata.obs['CellType'].str.startswith('fFB'), 'celltype'] = 'fFB'
adata.obs.loc[adata.obs['CellType'].str.startswith('dP'), 'celltype'] = 'dP'
adata.obs.loc[adata.obs['CellType'].str.startswith('dM'), 'celltype'] = 'dM'
adata.obs.loc[adata.obs['CellType'].str.startswith('dS'), 'celltype'] = 'dS'
adata.obs.loc[adata.obs['CellType'].str.startswith('Epi'), 'celltype'] = 'Epi'
adata.obs.loc[adata.obs['CellType'].str.startswith('DC'), 'celltype'] = 'DC'
adata.obs.loc[adata.obs['CellType'].str.startswith('Endo'), 'celltype'] = 'Endo'

ggs = ['CSF1', 'CSF2', 'XCL1', 'CCL5']

#reds = copy(matplotlib.cm.Reds)
#reds.set_under("lightgray")
fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(15, 10))
for i, pp in enumerate(['CellType'] + ggs + ['celltype']):
    ax = axs[i//3, i%3]
    ax.set_box_aspect(1)
    ss = 20 if i == 0 else 6
    alpha = 0.1 if i == 0 else 1
    ll = 'on data' if i < 5 else 'right margin'
    sc.pl.embedding(adata, basis='UMAP', color=pp, legend_loc=ll, legend_fontsize=10, show=False, ax=ax, size=ss, frameon=False, alpha=alpha)
fig.savefig('four_cytokines.pdf')
