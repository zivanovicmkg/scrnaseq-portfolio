#!/usr/bin/env python3
"""
Patch for Project 08: properly evaluate BBKNN.

The original notebook stored PCA under obsm['BBKNN'], so scib-metrics scored
the un-corrected PCA, not the BBKNN-modified graph. Fix:

1. Re-run BBKNN on the adata (modifies obsp directly)
2. Compute UMAP using the BBKNN-modified neighbors graph (now obsp has bbknn output)
3. Use diffusion-map embedding from the BBKNN graph as obsm['BBKNN'] so
   scib-metrics scores something that actually reflects the BBKNN integration.
   This is the approach used in the YosefLab scib-metrics lung tutorial.

Writes:
- updated hero_umap_4methods.png (BBKNN UMAP fixed)
- updated batch_composition_per_cluster.png (BBKNN cluster batch mix fixed)
- updated metrics_and_runtime.png
- updated integration_metrics_raw.csv + scaled.csv
- updated benchmark_summary.json
"""

import os, sys, time, json, warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse.csgraph import connected_components

warnings.filterwarnings('ignore')

sys.path.insert(0, '/home/marko-b2/upwork_portfolio/_shared/scripts')
try:
    from portfolio_utils import save_fig
    plt.style.use('/home/marko-b2/upwork_portfolio/_shared/styles/portfolio.mplstyle')
except Exception:
    def save_fig(fig, path, dpi=300):
        fig.savefig(path, dpi=dpi, bbox_inches='tight')

import bbknn
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

BASE = '/home/marko-b2/upwork_portfolio/08_integration_harmony'
FIG = f'{BASE}/figures'
RES = f'{BASE}/results'
NB_PATH = f'{BASE}/notebooks/analysis.ipynb'

BATCH_KEY = 'orig.ident'
LABEL_KEY = 'celltype_major'

# ---------- 1. Reload source + get existing Harmony/Scanorama embeddings ----------
# We need the embeddings we already computed. Easiest: reload Wu 2021 and re-run
# Harmony / Scanorama + add proper BBKNN. Since the notebook output is large
# and the embeddings are not saved separately, just re-run all three.
# This takes ~3 min on this dataset.

print('Loading Wu 2021 data...')
adata = sc.read('/home/marko-b2/upwork_portfolio/02_breast_cancer_TME/results/wu2021_3patient_TME_processed.h5ad')

# Strip prior outputs (same as notebook section 3)
for k in list(adata.obsm.keys()):
    if k.startswith('X_pca') or k.startswith('X_umap') or k in ['Unintegrated','Harmony','BBKNN','Scanorama']:
        del adata.obsm[k]
for u in ['neighbors', 'umap', 'pca', 'harmony', 'unintegrated', 'leiden']:
    if u in adata.uns:
        del adata.uns[u]
for p in list(adata.obsp.keys()):
    del adata.obsp[p]

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat', batch_key=BATCH_KEY)
print(f'HVGs: {int(adata.var.highly_variable.sum())}')

# ---------- 2. Unintegrated ----------
runtimes = {}
t0 = time.time()
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
adata.obsm['Unintegrated'] = adata.obsm['X_pca'].copy()
runtimes['Unintegrated'] = time.time() - t0

sc.pp.neighbors(adata, n_pcs=30, use_rep='Unintegrated', key_added='nbr_unint')
sc.tl.umap(adata, neighbors_key='nbr_unint')
adata.obsm['X_umap_unintegrated'] = adata.obsm['X_umap'].copy()
del adata.obsm['X_umap']

# ---------- 3. Harmony ----------
import harmonypy
t0 = time.time()
ho = harmonypy.run_harmony(adata.obsm['X_pca'].astype('float64'), adata.obs, BATCH_KEY, max_iter_harmony=20)
Z = ho.Z_corr
adata.obsm['Harmony'] = Z if Z.shape[0] == adata.n_obs else Z.T
runtimes['Harmony'] = time.time() - t0

sc.pp.neighbors(adata, n_pcs=30, use_rep='Harmony', key_added='nbr_harmony')
sc.tl.umap(adata, neighbors_key='nbr_harmony')
adata.obsm['X_umap_harmony'] = adata.obsm['X_umap'].copy()
del adata.obsm['X_umap']
print(f'Harmony: {runtimes["Harmony"]:.1f}s')

# ---------- 4. BBKNN (proper way) ----------
# bbknn writes its modified kNN graph to .obsp and .uns['neighbors']
ad_bbknn = adata.copy()
ad_bbknn.obsm['X_pca'] = adata.obsm['X_pca'].copy()

t0 = time.time()
bbknn.bbknn(ad_bbknn, batch_key=BATCH_KEY, n_pcs=30, neighbors_within_batch=10)
runtimes['BBKNN'] = time.time() - t0
print(f'BBKNN: {runtimes["BBKNN"]:.1f}s')

# UMAP from BBKNN graph (now ad_bbknn.obsp has bbknn output)
sc.tl.umap(ad_bbknn)
adata.obsm['X_umap_bbknn'] = ad_bbknn.obsm['X_umap'].copy()

# For scib-metrics Benchmarker we need an obsm embedding.
# Compute a diffusion-map embedding from the BBKNN graph - this propagates the
# batch-mixing through the connectivity into a continuous representation.
print('Computing diffmap from BBKNN graph for metric scoring...')
sc.tl.diffmap(ad_bbknn, n_comps=20)
# Drop the first DC (constant), keep 2..n
adata.obsm['BBKNN'] = ad_bbknn.obsm['X_diffmap'][:, 1:].copy()
print(f'BBKNN embedding (diffmap from bbknn graph): {adata.obsm["BBKNN"].shape}')

# Copy bbknn neighbors into ad for clustering later
adata.obsp['bbknn_connectivities'] = ad_bbknn.obsp['connectivities'].copy()
adata.obsp['bbknn_distances'] = ad_bbknn.obsp['distances'].copy()
adata.uns['bbknn_neighbors'] = {'connectivities_key': 'bbknn_connectivities',
                                 'distances_key': 'bbknn_distances',
                                 'params': {'n_neighbors': 30, 'method': 'umap'}}

# ---------- 5. Scanorama ----------
import scanorama
batches = sorted(adata.obs[BATCH_KEY].unique())
adatas_per_batch = [adata[adata.obs[BATCH_KEY] == b].copy() for b in batches]

t0 = time.time()
scanorama.integrate_scanpy(adatas_per_batch, dimred=50)
runtimes['Scanorama'] = time.time() - t0
print(f'Scanorama: {runtimes["Scanorama"]:.1f}s')

scanorama_embed = np.zeros((adata.n_obs, 50))
for ad_batch in adatas_per_batch:
    idx = adata.obs_names.get_indexer(ad_batch.obs_names)
    scanorama_embed[idx] = ad_batch.obsm['X_scanorama']
adata.obsm['Scanorama'] = scanorama_embed

sc.pp.neighbors(adata, n_pcs=30, use_rep='Scanorama', key_added='nbr_scanorama')
sc.tl.umap(adata, neighbors_key='nbr_scanorama')
adata.obsm['X_umap_scanorama'] = adata.obsm['X_umap'].copy()
del adata.obsm['X_umap']

# ---------- 6. Benchmark ----------
print('\nRunning Benchmarker on 4 embeddings...')
embedding_keys = ['Unintegrated', 'Harmony', 'BBKNN', 'Scanorama']

bm = Benchmarker(
    adata,
    batch_key=BATCH_KEY,
    label_key=LABEL_KEY,
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    embedding_obsm_keys=embedding_keys,
    n_jobs=4,
)
t0 = time.time()
bm.benchmark()
bench_time = time.time() - t0
print(f'Benchmark wall time: {bench_time:.1f}s')

results_df = bm.get_results(min_max_scale=False)
results_scaled = bm.get_results(min_max_scale=True)
results_df.to_csv(f'{RES}/integration_metrics_raw.csv')
results_scaled.to_csv(f'{RES}/integration_metrics_scaled.csv')
print('\nRaw metrics:')
print(results_df.round(3))

# scib-metrics ranking plot (overwrites previous)
for ext in ['svg', 'png']:
    p = f'{FIG}/scib_results.{ext}'
    if os.path.exists(p):
        os.remove(p)
    p2 = f'{FIG}/benchmark_ranking_table.{ext}'
    if os.path.exists(p2):
        os.remove(p2)
bm.plot_results_table(min_max_scale=True, show=False, save_dir=FIG)
for ext in ['svg', 'png']:
    candidate = f'{FIG}/scib_results.{ext}'
    if os.path.exists(candidate):
        os.rename(candidate, f'{FIG}/benchmark_ranking_table.{ext}')

# ---------- 7. Hero UMAP figure ----------
umap_keys = {
    'Unintegrated': 'X_umap_unintegrated',
    'Harmony':      'X_umap_harmony',
    'BBKNN':        'X_umap_bbknn',
    'Scanorama':    'X_umap_scanorama',
}

batch_vals = adata.obs[BATCH_KEY].astype(str)
label_vals = adata.obs[LABEL_KEY].astype(str)
unique_batches = sorted(batch_vals.unique())
unique_labels = sorted(label_vals.unique())
batch_pal = dict(zip(unique_batches, sns.color_palette('tab10', n_colors=len(unique_batches))))
label_pal = dict(zip(unique_labels, sns.color_palette('tab10', n_colors=len(unique_labels))))

fig, axes = plt.subplots(2, 4, figsize=(20, 10))
for col, (method, umap_key) in enumerate(umap_keys.items()):
    coords = adata.obsm[umap_key]
    ax = axes[0, col]
    for b in unique_batches:
        m = (batch_vals == b).values
        ax.scatter(coords[m, 0], coords[m, 1], s=2.5, c=[batch_pal[b]], label=b, alpha=0.6, linewidth=0)
    ax.set_title(f'{method}\n(colored by patient)', fontsize=11, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    if col == 3:
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=8, frameon=False, markerscale=3)
    ax = axes[1, col]
    for l in unique_labels:
        m = (label_vals == l).values
        ax.scatter(coords[m, 0], coords[m, 1], s=2.5, c=[label_pal[l]], label=l, alpha=0.6, linewidth=0)
    ax.set_title(f'{method}\n(colored by cell type)', fontsize=11, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    if col == 3:
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=8, frameon=False, markerscale=3)

plt.suptitle('Integration benchmark: UMAP comparison across 4 methods', fontsize=14, fontweight='bold', y=1.00)
plt.tight_layout()
save_fig(fig, f'{FIG}/hero_umap_4methods.png', dpi=300)
save_fig(fig, f'{FIG}/hero_umap_4methods.tiff', dpi=300)
plt.close()

# ---------- 8. Custom metric chart ----------
core_metrics = [c for c in ['KBET', 'iLISI', 'cLISI', 'Silhouette label', 'Silhouette batch', 'Graph connectivity'] if c in results_df.columns]
df_plot = results_df[core_metrics].copy()
if 'Metric Type' in df_plot.index:
    df_plot = df_plot.drop('Metric Type')
df_plot = df_plot.apply(pd.to_numeric, errors='coerce')

fig, axes = plt.subplots(1, 2, figsize=(15, 5), gridspec_kw={'width_ratios': [2.5, 1]})
method_palette = {'Unintegrated': '#bbbbbb', 'Harmony': '#1f77b4', 'BBKNN': '#ff7f0e', 'Scanorama': '#2ca02c'}
n_methods = len(df_plot.index)
n_metrics = len(core_metrics)
bar_w = 0.8 / n_methods
x = np.arange(n_metrics)
for i, method in enumerate(df_plot.index):
    vals = df_plot.loc[method, core_metrics].values.astype(float)
    axes[0].bar(x + i * bar_w, vals, width=bar_w, label=str(method),
                color=method_palette.get(str(method), 'gray'), edgecolor='white', linewidth=0.5)
axes[0].set_xticks(x + bar_w * (n_methods - 1) / 2)
axes[0].set_xticklabels(core_metrics, rotation=20, ha='right')
axes[0].set_ylabel('Metric value (higher = better)')
axes[0].set_title('Per-metric comparison', fontweight='bold')
axes[0].legend(loc='upper right', fontsize=9, frameon=False)
axes[0].grid(axis='y', alpha=0.3)

rt_methods = list(runtimes.keys())
rt_vals = [runtimes[m] for m in rt_methods]
axes[1].barh(rt_methods, rt_vals,
             color=[method_palette.get(m, 'gray') for m in rt_methods],
             edgecolor='white')
for i, (m, v) in enumerate(zip(rt_methods, rt_vals)):
    axes[1].text(v, i, f'  {v:.1f}s', va='center', fontsize=9)
axes[1].set_xlabel('Wall time (seconds, lower = faster)')
axes[1].set_title('Runtime', fontweight='bold')
axes[1].grid(axis='x', alpha=0.3)
plt.tight_layout()
save_fig(fig, f'{FIG}/metrics_and_runtime.png', dpi=300)
plt.close()

# ---------- 9. Per-cluster batch composition ----------
fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))
n_batches = len(unique_batches)
expected = 1.0 / n_batches

for col, method in enumerate(['Unintegrated', 'Harmony', 'BBKNN', 'Scanorama']):
    if method == 'BBKNN':
        # Cluster on BBKNN connectivities
        sc.tl.leiden(adata, resolution=0.5, neighbors_key='bbknn_neighbors',
                     obsp='bbknn_connectivities', key_added=f'leiden_{method}')
    else:
        nbr_key = {'Unintegrated': 'nbr_unint', 'Harmony': 'nbr_harmony', 'Scanorama': 'nbr_scanorama'}[method]
        sc.tl.leiden(adata, resolution=0.5, neighbors_key=nbr_key, key_added=f'leiden_{method}')
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

plt.suptitle('Per-cluster batch composition (uniform = well integrated)', fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
save_fig(fig, f'{FIG}/batch_composition_per_cluster.png', dpi=300)
plt.close()

# ---------- 10. Updated summary JSON ----------
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
    'bbknn_handling': 'BBKNN embedding for scib-metrics is the diffusion map computed from the bbknn-modified kNN graph; clustering for batch composition uses the bbknn graph directly.',
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
