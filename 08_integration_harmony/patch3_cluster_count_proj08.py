#!/usr/bin/env python3
"""
Mini-patch v3 for Project 08.

The BBKNN diffmap-based neighbors graph is dense within sub-populations, so
leiden at resolution=0.5 produces 35 clusters vs ~14 for other methods.
This makes the visual comparison unfair. Fix: tune BBKNN resolution down
to ~0.15 to produce a comparable cluster count.

This script only regenerates the batch composition plot - everything else
from patch2 stays as is.
"""

import os, sys, time, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')
sys.path.insert(0, '/home/marko-b2/upwork_portfolio/_shared/scripts')
try:
    from portfolio_utils import save_fig
    plt.style.use('/home/marko-b2/upwork_portfolio/_shared/styles/portfolio.mplstyle')
except Exception:
    def save_fig(fig, path, dpi=300):
        fig.savefig(path, dpi=dpi, bbox_inches='tight')

import harmonypy
import bbknn
import scanorama

BASE = '/home/marko-b2/upwork_portfolio/08_integration_harmony'
FIG = f'{BASE}/figures'
RES = f'{BASE}/results'

BATCH_KEY = 'orig.ident'
LABEL_KEY = 'celltype_major'

print('Rebuilding integrations for batch composition plot...')
adata = sc.read('/home/marko-b2/upwork_portfolio/02_breast_cancer_TME/results/wu2021_3patient_TME_processed.h5ad')

for k in list(adata.obsm.keys()):
    if k.startswith('X_') or k in ['Unintegrated','Harmony','BBKNN','Scanorama']:
        del adata.obsm[k]
for u in ['neighbors','umap','pca','harmony','unintegrated','leiden']:
    if u in adata.uns:
        del adata.uns[u]
for p in list(adata.obsp.keys()):
    del adata.obsp[p]

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat', batch_key=BATCH_KEY)
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
adata.obsm['Unintegrated'] = adata.obsm['X_pca'].copy()
sc.pp.neighbors(adata, n_pcs=30, use_rep='Unintegrated', key_added='nbr_unint')

ho = harmonypy.run_harmony(adata.obsm['X_pca'].astype('float64'), adata.obs, BATCH_KEY, max_iter_harmony=20)
Z = ho.Z_corr
adata.obsm['Harmony'] = Z if Z.shape[0] == adata.n_obs else Z.T
sc.pp.neighbors(adata, n_pcs=30, use_rep='Harmony', key_added='nbr_harmony')

ad_bbknn = adata.copy()
ad_bbknn.obsm['X_pca'] = adata.obsm['X_pca'].copy()
bbknn.bbknn(ad_bbknn, batch_key=BATCH_KEY, n_pcs=30, neighbors_within_batch=10)
sc.tl.diffmap(ad_bbknn, n_comps=20)
adata.obsm['BBKNN'] = ad_bbknn.obsm['X_diffmap'][:, 1:].copy()
sc.pp.neighbors(adata, n_pcs=adata.obsm['BBKNN'].shape[1], use_rep='BBKNN', key_added='nbr_bbknn')

batches = sorted(adata.obs[BATCH_KEY].unique())
adatas_per_batch = [adata[adata.obs[BATCH_KEY] == b].copy() for b in batches]
scanorama.integrate_scanpy(adatas_per_batch, dimred=50)
scanorama_embed = np.zeros((adata.n_obs, 50))
for ab in adatas_per_batch:
    idx = adata.obs_names.get_indexer(ab.obs_names)
    scanorama_embed[idx] = ab.obsm['X_scanorama']
adata.obsm['Scanorama'] = scanorama_embed
sc.pp.neighbors(adata, n_pcs=30, use_rep='Scanorama', key_added='nbr_scanorama')

# ---------- Cluster with tuned resolution per method ----------
# Goal: ~12-15 clusters for each method.
batch_vals = adata.obs[BATCH_KEY].astype(str)
unique_batches = sorted(batch_vals.unique())
batch_pal = dict(zip(unique_batches, sns.color_palette('tab10', n_colors=len(unique_batches))))
expected = 1.0 / len(unique_batches)

method_config = {
    'Unintegrated': ('nbr_unint',     0.5),
    'Harmony':      ('nbr_harmony',   0.5),
    'BBKNN':        ('nbr_bbknn',     0.15),  # tuned down: diffmap-based graph is dense
    'Scanorama':    ('nbr_scanorama', 0.5),
}

print('\nLeiden clustering per method:')
for method, (nbr_key, res) in method_config.items():
    sc.tl.leiden(adata, resolution=res, neighbors_key=nbr_key, key_added=f'leiden_{method}')
    n_clust = adata.obs[f'leiden_{method}'].nunique()
    print(f'  {method:<15} resolution={res:<5} -> {n_clust} clusters')

# ---------- Plot ----------
fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))
for col, method in enumerate(['Unintegrated', 'Harmony', 'BBKNN', 'Scanorama']):
    clusters = adata.obs[f'leiden_{method}'].astype(str)
    df_obs = pd.DataFrame({'cluster': clusters.values, 'batch': batch_vals.values})
    comp = df_obs.groupby('cluster')['batch'].value_counts(normalize=True).unstack(fill_value=0)
    comp = comp.reindex(columns=unique_batches, fill_value=0)
    comp = comp.sort_index(key=lambda x: x.astype(int) if x.str.isnumeric().all() else x)
    
    ax = axes[col]
    bottom = np.zeros(len(comp))
    for b in unique_batches:
        ax.bar(range(len(comp)), comp[b].values, bottom=bottom,
               color=batch_pal[b], edgecolor='white', linewidth=0.5, label=b)
        bottom = bottom + comp[b].values
    ax.axhline(expected, color='red', linestyle='--', lw=1, alpha=0.7, label=f'expected ({expected:.2f})')
    ax.set_xlabel(f'Leiden cluster (n={len(comp)})')
    ax.set_ylabel('Batch fraction' if col == 0 else '')
    ax.set_title(method, fontweight='bold')
    ax.set_ylim(0, 1)
    if col == 3:
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=8, frameon=False)

plt.suptitle('Per-cluster batch composition (uniform = well integrated)',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
save_fig(fig, f'{FIG}/batch_composition_per_cluster.png', dpi=300)
plt.close()
print(f'\nSaved: {FIG}/batch_composition_per_cluster.png')
print('DONE.')
