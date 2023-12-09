import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors, cm
plt.rcParams["pdf.use14corefonts"] = True

rpkm = pd.read_csv('table_all_genes_paired_comparisons.csv', index_col=0)
gty = pd.read_csv('../genomes_annotations/gencode.v37.geneid.genename.genetype.txt', sep='\t', header=None, index_col=0)
gty.index = gty.index.str.replace('\\..*', '')

##S1sp-specific upregulated genes
ind1 = ((rpkm['DE_summary'].str.count('Up-reg') == 1) & rpkm['DE_summary'].str.contains('S1sp:Up-reg'))
ind2 = (gty.loc[rpkm['Ensembl ID'], 1] == 'protein_coding')
ggs = rpkm.index[ind1.values & ind2.values]

##calculate log2-transformed fold change between C2+ and C1+HLA-C for each donor
ind = rpkm.columns.str.startswith('Y')
rpkm = rpkm.loc[:, ind]
gps = ['S1sp', 'L1sp', 'S1L1dp']
mats = []
for gp in gps:
    mat = rpkm.loc[:, rpkm.columns.str.contains(gp)]
    mat = mat.iloc[:, 3:].to_numpy()/mat.iloc[:, :3].to_numpy()
    mats.append(mat)
mat = np.concatenate(mats, axis=1)
mat[mat > 10] = 10
mat[np.isnan(mat)] = 1
mat = np.log2(mat)
mat[mat < -np.log2(10)] = -np.log2(10)
mat = pd.DataFrame(mat, index=rpkm.index)

mat = mat.loc[ggs]
vmin = np.percentile(mat, 5)
vmax = np.percentile(mat, 95)
lgg = ['CSF2', 'XCL1', 'CCL1', 'CCL3L3', 'CCL4', 'CCL4L2', 'CCL3', 'CXCL13', 'XCL2', 'IL3', 'IL31']

cols = list(map(colors.to_hex, cm.tab10.colors))
pcols = [cols[1]]*3 + [cols[0]]*3 + [cols[2]]*3
ax = sns.clustermap(mat, row_cluster=True, col_cluster=False, col_colors=pcols, vmin=vmin, vmax=vmax, center=0, yticklabels=False, xticklabels=False, figsize=(3, 8), cmap='RdBu_r', colors_ratio=0.02, method='ward', metric='euclidean')
ind = ax.dendrogram_row.reordered_ind
ax.ax_heatmap.set_yticks([np.where(mat.index[ind]==g)[0][0] + 0.5 for g in lgg])
ax.ax_heatmap.set_yticklabels(lgg)
ax.ax_heatmap.set_xticks([1, 4, 7])
ax.ax_heatmap.set_xticklabels(gps)
ax.ax_heatmap.tick_params(axis='both', labelsize=5)
ax.savefig('S1sp_up_genes_lgfd_heatmap.pdf')
