
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc



sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_versions()
#results_file = './write/pa.h5ad'
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(3, 3), facecolor='white')

raw_UMIcounts = pd.read_table("./data/GSM4798908_B2019-1.expression_matrix.txt", header=0, index_col=0)
#raw_UMIcounts.head(10)

raw_UMIcounts_pl = pd.read_table("./data/GSM4798909_B2019-2.expression_matrix.txt", header=0, index_col=0)

raw_UMIcounts_nl= pd.read_table("./data/GSM4798910_B2019-3.expression_matrix.txt", header=0, index_col=0)

pi_gene_bar= raw_UMIcounts.columns

pi_gene_bar=pd.Series(pi_gene_bar)+"_primary"

raw_UMIcounts.columns= pi_gene_bar

raw_UMIcounts_pl.columns = pd.Series(raw_UMIcounts_pl.columns)+"_positive"

raw_UMIcounts_nl.columns = pd.Series(raw_UMIcounts_nl.columns)+"_negtive"

counts_join=raw_UMIcounts.join(raw_UMIcounts_pl,how="inner")

counts_join= counts_join.join(raw_UMIcounts_nl,how="inner")



counts_join.T.to_csv("./data/counts_join.csv")


adata = sc.read_csv("./data/counts_join.csv", first_column_names=True)


adata.write('./data/join_adata.h5ad')
