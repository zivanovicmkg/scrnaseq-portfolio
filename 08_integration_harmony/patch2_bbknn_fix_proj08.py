#!/usr/bin/env python3
"""
Mini-patch for Project 08 patch script.

The main patch ran successfully through the Benchmarker (Total scores OK).
It crashed at the per-cluster batch composition step due to a scanpy API:
sc.tl.leiden does not accept both neighbors_key and obsp at the same time.

This script:
1. Re-runs all 4 integrations (fast: ~3 min total) to rebuild adata in memory
2. Generates the batch composition plot with the correct leiden call
3. Saves benchmark_summary.json (which the previous patch did not reach)

The hero UMAP and metrics chart were already saved by the previous patch.
"""

import os, sys, time, json, warnings
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
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

BASE = '/home/marko-b2/upwork_portfolio/08_integration_harmony'
FIG = f'{BASE}/figures'
RES = f'{BASE}/results'

BATCH_KEY = 'orig.ident'
LABEL_KEY = 'celltype_major'

print('Loading + rebuilding all integrations...')
adata = sc.read('/home/marko-b2/upwork_portfolio/02_breast_cancer_TME/results/wu2021_3patient_TME_processed.h5ad')

# Strip prior
for k in list(adata.obsm.keys()):
    if k.startswith('X_') or k in ['Unintegrated','Harmony','BBKNN','Scanorama']:
        del adata.obsm[k]
for u in ['neighbors','umap','pca','harmony','unintegrated','leiden']:
    if u in adata.uns:
        del adata.uns[u]
for p in list(adata.obsp.keys()):
    del adata.obsp[p]

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat', batch_key=BATCH_KEY)

runtimes = {}
# Unintegrated
t0 = time.time()
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
adata.obsm['Unintegrated'] = adata.obsm['X_pca'].copy()
runtimes['Unintegrated'] = time.time() - t0
sc.pp.neighbors(adata, n_pcs=30, use_rep='Unintegrated', key_added='nbr_unint')

# Harmony
t0 = time.time()
ho = harmonypy.run_harmony(adata.obsm['X_pca'].astype('float64'), adata.obs, BATCH_KEY, max_iter_harmony=20)
Z = ho.Z_corr
adata.obsm['Harmony'] = Z if Z.shape[0] == adata.n_obs else Z.T
runtimes['Harmony'] = time.time() - t0
sc.pp.neighbors(adata, n_pcs=30, use_rep='Harmony', key_added='nbr_harmony')

# BBKNN
ad_bbknn = adata.copy()
ad_bbknn.obsm['X_pca'] = adata.obsm['X_pca'].copy()
t0 = time.time()
bbknn.bbknn(ad_bbknn, batch_key=BATCH_KEY, n_pcs=30, neighbors_within_batch=10)
runtimes['BBKNN'] = time.time() - t0
sc.tl.diffmap(ad_bbknn, n_comps=20)
adata.obsm['BBKNN'] = ad_bbknn.obsm['X_diffmap'][:, 1:].copy()
# Re-run neighbors ON the BBKNN diffmap so leiden has a proper nbr key (no obsp conflict)
sc.pp.neighbors(adata, n_pcs=adata.obsm['BBKNN'].shape[1], use_rep='BBKNN', key_added='nbr_bbknn')

# Scanorama
batches = sorted(adata.obs[BATCH_KEY].unique())
adatas_per_batch = [adata[adata.obs[BATCH_KEY] == b].copy() for b in batches]
t0 = time.time()
scanorama.integrate_scanpy(adatas_per_batch, dimred=50)
runtimes['Scanorama'] = time.time() - t0
scanorama_embed = np.zeros((adata.n_obs, 50))
for ab in adatas_per_batch:
    idx = adata.obs_names.get_indexer(ab.obs_names)
    scanorama_embed[idx] = ab.obsm['X_scanorama']
adata.obsm['Scanorama'] = scanorama_embed
sc.pp.neighbors(adata, n_pcs=30, use_rep='Scanorama', key_added='nbr_scanorama')

print(f'\nRuntimes: {runtimes}')

# ---------- Benchmark (we need results_df again for the summary JSON) ----------
embedding_keys = ['Unintegrated', 'Harmony', 'BBKNN', 'Scanorama']
bm = Benchmarker(adata, batch_key=BATCH_KEY, label_key=LABEL_KEY,
                 bio_conservation_metrics=BioConservation(),
                 batch_correction_metrics=BatchCorrection(),
                 embedding_obsm_keys=embedding_keys, n_jobs=4)
t0 = time.time()
bm.benchmark()
bench_time = time.time() - t0
results_df = bm.get_results(min_max_scale=False)
print(f'\nBenchmark done in {bench_time:.1f}s')

# ---------- Per-cluster batch composition (FIXED: only neighbors_key, no obsp) ----------
print('\nGenerating per-cluster batch composition plot...')
batch_vals = adata.obs[BATCH_KEY].astype(str)
unique_batches = sorted(batch_vals.unique())
batch_pal = dict(zip(unique_batches, sns.color_palette('tab10', n_colors=len(unique_batches))))

fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))
expected = 1.0 / len(unique_batches)
nbr_map = {
    'Unintegrated': 'nbr_unint',
    'Harmony':      'nbr_harmony',
    'BBKNN':        'nbr_bbknn',
    'Scanorama':    'nbr_scanorama',
}

for col, method in enumerate(['Unintegrated', 'Harmony', 'BBKNN', 'Scanorama']):
    sc.tl.leiden(adata, resolution=0.5,
                 neighbors_key=nbr_map[method],
                 key_added=f'leiden_{method}')
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
    ax.set_xlabel('Leiden cluster')
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
print(f'Saved: {FIG}/batch_composition_per_cluster.png')

# ---------- Save summary JSON ----------
summary = {
    'dataset': 'Wu 2021 breast cancer, 3 patients, 12,962 cells',
    'batch_key': BATCH_KEY,
    'label_key': LABEL_KEY,
    'n_batches': int(adata.obs[BATCH_KEY].nunique()),
    'n_labels': int(adata.obs[LABEL_KEY].nunique()),
    'n_hvg': int(adata.var.highly_variable.sum()),
    'methods': embedding_keys,
    'runtimes_sec': {k: float(v) for k, v in runtimes.items()},
    'benchmark_time_sec': float(bench_time),
    'bbknn_handling': 'BBKNN embedding for scib-metrics is the diffusion map computed from the bbknn-modified kNN graph.',
}
for agg in ['Total', 'Batch correction', 'Bio conservation']:
    if agg in results_df.columns:
        summary[f'aggregate_{agg.lower().replace(" ", "_")}'] = {
            str(m): float(pd.to_numeric(results_df.loc[m, agg], errors='coerce'))
            for m in results_df.index if m != 'Metric Type'
        }
with open(f'{RES}/benchmark_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print('\n=== AGGREGATE SCORES ===')
for agg in ['Total', 'Batch correction', 'Bio conservation']:
    if agg in results_df.columns:
        print(f'\n{agg}:')
        for m in results_df.index:
            if m != 'Metric Type':
                v = pd.to_numeric(results_df.loc[m, agg], errors='coerce')
                print(f'  {str(m):<15} {v:.3f}')
print('\nDONE.')
