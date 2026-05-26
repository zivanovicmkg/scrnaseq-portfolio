#!/usr/bin/env python3
"""
Patch + rerun comparison: fix marker Jaccard NaN issue.
Issue: py_markers['cluster'] is string, r_markers['cluster'] is integer.
The query() calls with different types returned empty sets, so all Jaccard = NaN.
"""
import json
from pathlib import Path

NB_PATH = '/home/marko-b2/upwork_portfolio/10_crossplatform_reanalysis/notebooks/03_comparison.ipynb'

with open(NB_PATH) as f:
    nb = json.load(f)

# Replace the marker-Jaccard cell with a type-safe version
new_marker_src = """# For each Python cluster, find the R cluster with maximum cell overlap
py_to_r_match = {}
for py_cl in ct_mat.index:
    best_r_cl = ct_mat.loc[py_cl].idxmax()
    py_to_r_match[py_cl] = best_r_cl

print('Python cluster -> R cluster (by max cell overlap):')
for k, v in py_to_r_match.items():
    print(f'  {k} -> {v}  (n={ct_mat.loc[k, v]})')

# Normalize both marker DataFrames to string cluster column for safe matching
py_markers['cluster'] = py_markers['cluster'].astype(str)
r_markers['cluster']  = r_markers['cluster'].astype(str)

# Jaccard at top 10, 25, 50
jaccard_rows = []
for n_top in [10, 25, 50]:
    for py_cl, r_cl in py_to_r_match.items():
        py_cl_str = str(py_cl)
        r_cl_str = str(r_cl)
        py_genes = set(py_markers[py_markers['cluster'] == py_cl_str].head(n_top)['gene'])
        r_genes  = set(r_markers[r_markers['cluster'] == r_cl_str].head(n_top)['gene'])
        if not py_genes or not r_genes:
            jacc = float('nan')
        else:
            jacc = len(py_genes & r_genes) / len(py_genes | r_genes)
        jaccard_rows.append({'top_n': n_top, 'py_cluster': py_cl, 'r_cluster': r_cl,
                              'n_py_genes': len(py_genes), 'n_r_genes': len(r_genes),
                              'jaccard': jacc})
jacc_df = pd.DataFrame(jaccard_rows)
jacc_df.to_csv(f'{RES_DIR}/marker_jaccard.csv', index=False)

print('\\nMean Jaccard across matched clusters:')
for n_top in [10, 25, 50]:
    mean_j = jacc_df.query(f'top_n == {n_top}')['jaccard'].mean()
    print(f'  top {n_top}: {mean_j:.3f}')

print('\\nPer-cluster Jaccard (top 10):')
print(jacc_df.query('top_n == 10')[['py_cluster','r_cluster','n_py_genes','n_r_genes','jaccard']])
"""

# Find and patch the marker-Jaccard cell
patched_jacc = False
for c in nb['cells']:
    if c['cell_type'] == 'code':
        src = ''.join(c['source']) if isinstance(c['source'], list) else c['source']
        if 'py_to_r_match' in src and 'jaccard' in src.lower() and 'jacc_df' in src:
            c['source'] = new_marker_src.splitlines(keepends=True)
            c['outputs'] = []
            c['execution_count'] = None
            patched_jacc = True
            print('Patched marker-Jaccard cell')
            break

# Also patch the comparison_summary cell to handle NaN gracefully
new_summary_src = """def safe_float(x):
    try:
        v = float(x)
        return v if v == v else None  # NaN check (NaN != NaN)
    except (ValueError, TypeError):
        return None

comparison_summary = {
    'dataset': 'PBMC 3k (10X Genomics)',
    'cells_common': int(len(common)),
    'cells_python': int(py_sum['n_cells_after_qc']),
    'cells_r': int(r_sum['n_cells_after_qc']),
    'clusters_python': int(py_sum['n_leiden_clusters']),
    'clusters_r': int(r_sum['n_seurat_clusters']),
    'cluster_agreement': {
        'ARI_leiden_vs_seurat': float(ari),
        'NMI_leiden_vs_seurat': float(nmi),
        'ARI_celltype': float(ari_ct),
        'NMI_celltype': float(nmi_ct),
    },
    'pca_agreement': {
        'mean_abs_spearman_top10': float(matched_corrs[:10].mean()),
        'mean_abs_spearman_all': float(matched_corrs.mean()),
        'per_pc': [{'py_pc': int(p+1), 'r_pc': int(r+1), 'abs_rho': float(c)}
                   for p, r, c in zip(row_ind, col_ind, matched_corrs)],
    },
    'umap_agreement': {
        'procrustes_disparity': float(disp),
    },
    'marker_agreement': {
        'mean_jaccard_top10': safe_float(jacc_df.query('top_n == 10').jaccard.mean()),
        'mean_jaccard_top25': safe_float(jacc_df.query('top_n == 25').jaccard.mean()),
        'mean_jaccard_top50': safe_float(jacc_df.query('top_n == 50').jaccard.mean()),
    },
    'runtime_sec': {
        'scanpy_total': float(py_sum['runtime_sec']['total']),
        'seurat_total': float(r_sum['runtime_sec']['total']),
    },
}

with open(f'{RES_DIR}/comparison_summary.json', 'w') as f:
    json.dump(comparison_summary, f, indent=2, default=str)

print(json.dumps({k: v for k, v in comparison_summary.items() if k != 'pca_agreement'}, indent=2, default=str))
print('\\nSaved:', f'{RES_DIR}/comparison_summary.json')
"""

patched_sum = False
for c in nb['cells']:
    if c['cell_type'] == 'code':
        src = ''.join(c['source']) if isinstance(c['source'], list) else c['source']
        if 'comparison_summary' in src and 'mean_jaccard_top10' in src:
            c['source'] = new_summary_src.splitlines(keepends=True)
            c['outputs'] = []
            c['execution_count'] = None
            patched_sum = True
            print('Patched summary cell')
            break

# Clear outputs from all cells onward to re-run cleanly
for c in nb['cells']:
    if c['cell_type'] == 'code':
        c['outputs'] = []
        c['execution_count'] = None

with open(NB_PATH, 'w') as f:
    json.dump(nb, f, indent=1)

print(f'\\nDone. Patched marker cell: {patched_jacc}, summary cell: {patched_sum}')
