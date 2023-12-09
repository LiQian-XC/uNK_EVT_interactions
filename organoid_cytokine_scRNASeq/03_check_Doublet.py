import scanpy as sc
import scrublet as scr
import pandas as pd
import numpy as np
import sys
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

dd = '/home/ql312/rds/rds-turco-lab2-jxGoj1xLQV4/trophoblast_organoid_from_Sanger/'
sp = sys.argv[1]
print(sp)

##load data
adata = sc.read_10x_mtx(f'{dd}/{sp}/cellranger_res/{sp}/outs/filtered_feature_bc_matrix/')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
info1 = pd.read_csv(f'{dd}/{sp}/souporcell_res/clusters.tsv', delimiter='\t')
adata.obs['status'] = info1['status'].values
adata.obs['assignment'] = info1['assignment'].values
info2 = pd.read_csv('Troph_Org_metadata.csv')
adata.obs[info2.columns] = info2[info2['sample']==sp].iloc[0].values
cq = '1' if sp.startswith('6044') else '2'
info3 = pd.read_csv(f'souporcell_shared_samples_{cq}.csv', index_col=0)
adata.obs['ID'] = 'unassigned'
for i in range(3):
	ind = (adata.obs['status']=='singlet') & (adata.obs['assignment']==str(i))
	adata.obs.loc[ind, 'ID'] = info3.index[info3[sp]==i][0]

##doublet score from scrublet
scrub = scr.Scrublet(adata.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
adata.obs['doublet_scores'] = doublet_scores
adata.obs['predicted_doublets'] = predicted_doublets

##cluster
sc.pp.filter_genes(adata, min_cells=20)
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, layer='counts', n_top_genes=2000, flavor='seurat_v3', subset=True)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=40)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata)
adata.obs['fleiden'] = adata.obs['leiden']

##subcluster
cls = np.unique(adata.obs['leiden'])
for cl in cls:
	ss = sum(adata.obs['leiden'] == cl)
	rr = np.maximum(np.maximum(np.log10(ss)-1, 0)**2, 0.1)
	sc.tl.leiden(adata, restrict_to=('leiden', [cl]), key_added='leiden_sub', resolution=rr)
	adata.obs['leiden'] = adata.obs['leiden_sub']

##calculate significance
def pval_fun(vals):
	med = np.median(vals)
	mad = np.median(vals[vals > med] - med)*1.4826
	tvals = (vals - med)/mad
	pvals = 1 - norm.cdf(tvals, loc=0, scale=1)
	bf_pvals = multipletests(pvals, method='fdr_bh')[1]
	return bf_pvals

med_vals = []
for cl in np.unique(adata.obs['leiden']):
	ind = adata.obs['leiden']==cl
	vv = np.median(adata.obs.loc[ind, 'doublet_scores'])
	med_vals.append(vv)
	adata.obs.loc[ind, 'median_doublet_scores'] = vv
med_vals = np.array(med_vals)
bf_pvals = pval_fun(med_vals)
for i, cl in enumerate(np.unique(adata.obs['leiden'])):
	ind = adata.obs['leiden'] == cl
	adata.obs.loc[ind, 'doublet_bf_pvals'] = bf_pvals[i]
adata.obs['doublet_local_pred'] = adata.obs['doublet_bf_pvals'] < 0.05
adata.obs['doublet_local_pred'] = adata.obs['doublet_local_pred'].astype(str)
adata.obs.to_csv(f'sample_info_{sp}.csv')
