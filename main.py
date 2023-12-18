import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc


adata=sc.read_h5ad("./data/join_adata.h5ad")
sc.pl.highest_expr_genes(adata, n_top=20, )

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=5)





adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

adata = adata[adata.obs.n_genes_by_counts <3000, :]
adata = adata[adata.obs.pct_counts_mt <20, :]

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

a=[]
for i in adata.obs.index:
    a.append(i[13:])

from collections import Counter



count = Counter(a)
print(count)

adata.obs["disease_group"]=a


sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)



sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata)

adata = adata[:, adata.var.highly_variable]


Counter(adata.var.highly_variable)



sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True)


results_file= './data/first_ana.h5ad'
adata.write(results_file)

adata=sc.read_h5ad("./data/first_ana.h5ad")


sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)

sc.tl.umap(adata,min_dist=0.3)
sc.pl.umap(adata, color="disease_group",)






from matplotlib.pyplot import rc_context

with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color='disease_group', add_outline=True,
               legend_fontsize=12, legend_fontoutline=2,frameon=False,
               title='clustering of cells', palette='Paired')

sc.tl.leiden(adata)

sc.tl.louvain(adata,resolution=0.6)



sc.pl.umap(adata, color=['leiden','louvain'], add_outline=True,legend_loc='on data',
            legend_fontsize=12, legend_fontoutline=2,frameon=False,
            title=['clustering of cells','clustering of cells'], palette='Set1')


# sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
#
#
# result = adata.uns['rank_genes_groups']
# groups = result['names'].dtype.names
# marker=pd.DataFrame(
#     {group + '_' + key[:1]: result[key][group]
#     for group in groups for key in ['names', 'pvals']})
#
# marker.to_csv("./data/marker_louvain.csv")

cluster2annotation = {
    '0': 'Bcells',
    '1': 'myofibroblasts',
    '2': 'NaiveT',
    '3': 'NaiveT',
    '4': 'cancercells',
    '5': 'Bcells',
    '6': 'NaiveT',
    '7': 'NaiveT',
    '8': 'myofibroblasts',
    '9': 'CD8Effector',
    '10': 'myeloid',
    '11': 'Bcells',
    '12': 'endothelial',
    '13': 'Macrophages',
    '14': 'cancerstemcells',
    '15': 'myeloid',
    '16': 'CXCL14cancer',
    '17': 'plasma',
    '18': 'matureDC',
    '19': 'pericytes'

}

# add a new `.obs` column called `cell type` by mapping clusters to annotation using pandas `map` function
adata.obs['cell type'] = adata.obs['louvain'].map(cluster2annotation).astype('category')

marker_genes_dict = {
    'cancercells': ['KRT19'],
    'cancerstemcells': ['KRT19','TOP2A'],
    'CXCL14cancer': ['CXCL14'],
    #'NaiveT': [],
    #'CD8Effector': [],
    'Bcells': ['CD79A','CD79B'],
    'Macrophages': ['LYZ','IL1B'],
    'myeloid':['LYZ'],
    'matureDC':['LAMP3','CCR7'],
    'plasma':['JCHAIN','IGHG3','MZB1'],
    'endothelial':['MCAM','PECAM1'],
    'pericytes':['ACTA2','TAGLN','MCAM'],
    'myofibroblasts':['LUM','DCN','TAGLN'],
    'cyclingcells':['TOP2A']
}



sc.pl.dotplot(adata, marker_genes_dict, 'louvain', dendrogram=True)


ax = sc.pl.stacked_violin(adata, marker_genes_dict, groupby='louvain', swap_axes=False, dendrogram=True,cmap='Paired_r')

sc.pl.tracksplot(adata, marker_genes_dict, groupby='louvain', swap_axes=False, dendrogram=True,cmap='Paired_r')



with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(adata, ['CD79A', 'CD79B'], groupby='louvain' )


sc.pl.umap(adata, color='louvain', legend_loc='on data',
           frameon=False, legend_fontsize=10, legend_fontoutline=2,title="")


sc.pl.umap(adata, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=5, legend_fontoutline=0.5,save="jmzeng.jpg")

adata.obs['cell type'].cat.categories


adata_new=adata[(adata.obs['cell type']== 'CXCL14cancer')|(adata.obs['cell type']=='cancercells')|(adata.obs['cell type']== 'cancerstemcells'), :]

adata_new.obs['cell type'].cat.categories



sc.tl.paga(adata_new, groups='cell type')

sc.pl.paga(adata_new, color=['cell type','CXCL14'])

sc.tl.draw_graph(adata_new, init_pos='paga')


adata_new.uns['iroot'] = np.flatnonzero(adata_new.obs['cell type']  == 'cancerstemcells')[0]
sc.tl.dpt(adata_new)



sc.pl.draw_graph(adata_new, color=['cell type', 'dpt_pseudotime','disease_group'], legend_loc='on data')


sc.tl.paga(adata, groups='cell type')
sc.pl.paga(adata, color=['cell type'],node_size_scale=0.5)

sc.tl.draw_graph(adata, init_pos='paga')


sc.pl.draw_graph(adata, color=['cell type'], legend_loc='on data',legend_fontsize='xx-small',legend_fontweight='normal')

adata.uns['iroot'] = np.flatnonzero(adata.obs['cell type']  == 'NaiveT')[0]
sc.tl.dpt(adata)


sc.pl.draw_graph(adata, color=['cell type', 'dpt_pseudotime','disease_group'], legend_loc='on data',legend_fontsize='xx-small',legend_fontweight='normal')
















