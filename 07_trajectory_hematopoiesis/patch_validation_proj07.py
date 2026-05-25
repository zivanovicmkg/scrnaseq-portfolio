#!/usr/bin/env python3
"""
Patch script for Project 07: re-compute per-branch validation with label mapping.

Mapping (our label -> published label(s)):
  Ery   -> Ery
  Mono  -> Mono
  DC    -> max(pDC, cDC) per cell  (sum of probs)
  Lymph -> CLP

Reads existing per-cell results + raw h5ad, writes updated summary JSON + plot.
Does NOT re-run Palantir.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import warnings

warnings.filterwarnings('ignore')

# Portfolio style
sys.path.insert(0, '/home/marko-b2/upwork_portfolio/_shared/scripts')
try:
    plt.style.use('/home/marko-b2/upwork_portfolio/_shared/styles/portfolio.mplstyle')
except Exception:
    pass

BASE = '/home/marko-b2/upwork_portfolio/07_trajectory_hematopoiesis'
DATA = f'{BASE}/data'
FIG = f'{BASE}/figures'
RES = f'{BASE}/results'

# ---------- 1. Load published data ----------
print('Loading published h5ad...')
ad = sc.read(f'{DATA}/human_cd34_bm_rep1.h5ad')
orig_pt = ad.obs['palantir_pseudotime']
orig_dp = ad.obs['palantir_diff_potential']
orig_branch = ad.obsm['palantir_branch_probs']
orig_ct = list(ad.uns['palantir_branch_probs_cell_types'])
orig_ct = [str(x) for x in orig_ct]  # strip np.str_
print(f'Published terminals: {orig_ct}')

orig_branch_df = pd.DataFrame(orig_branch, index=ad.obs_names, columns=orig_ct)

# ---------- 2. Load our per-cell results ----------
print('Loading our per-cell results...')
our = pd.read_csv(f'{RES}/palantir_per_cell.csv', index_col=0)
our_terms = [c.replace('fate_', '') for c in our.columns if c.startswith('fate_')]
print(f'Our terminals: {our_terms}')

# Align cells
common = our.index.intersection(orig_branch_df.index)
our_aligned = our.loc[common]
orig_pt_aligned = orig_pt.loc[common]
orig_dp_aligned = orig_dp.loc[common]
orig_branch_aligned = orig_branch_df.loc[common]
print(f'Common cells: {len(common)}')

# ---------- 3. Label mapping ----------
# Build per-cell published probability for each of our labels
published_for_our_label = {}

# Ery -> Ery (direct)
if 'Ery' in orig_branch_aligned.columns:
    published_for_our_label['Ery'] = orig_branch_aligned['Ery'].values

# Mono -> Mono (direct)
if 'Mono' in orig_branch_aligned.columns:
    published_for_our_label['Mono'] = orig_branch_aligned['Mono'].values

# DC -> pDC + cDC (sum, since both are DC subtypes lumped in our analysis)
dc_components = [c for c in ['pDC', 'cDC'] if c in orig_branch_aligned.columns]
if dc_components:
    published_for_our_label['DC'] = orig_branch_aligned[dc_components].sum(axis=1).values
    print(f'DC mapped to sum of: {dc_components}')

# Lymph -> CLP (direct)
if 'CLP' in orig_branch_aligned.columns:
    published_for_our_label['Lymph'] = orig_branch_aligned['CLP'].values
    print('Lymph mapped to: CLP')

# ---------- 4. Compute Spearman per branch ----------
branch_rhos = {}
for our_label in our_terms:
    if our_label not in published_for_our_label:
        print(f'WARN: no published mapping for {our_label}')
        continue
    our_vals = our_aligned[f'fate_{our_label}'].values
    pub_vals = published_for_our_label[our_label]
    rho, pval = spearmanr(our_vals, pub_vals)
    branch_rhos[our_label] = {
        'spearman_rho': float(rho),
        'p': float(pval),
        'mapped_to_published': dc_components if our_label == 'DC' else (
            'CLP' if our_label == 'Lymph' else our_label
        ),
    }
    map_str = branch_rhos[our_label]['mapped_to_published']
    print(f'  {our_label:<6} (-> {map_str})  rho = {rho:.3f}  (p = {pval:.2e})')

# ---------- 5. Pseudotime + entropy vs published ----------
rho_pt, p_pt = spearmanr(our_aligned['palantir_pseudotime'], orig_pt_aligned)
rho_dp, p_dp = spearmanr(our_aligned['palantir_entropy'], orig_dp_aligned)
print(f'\nPseudotime vs published:       rho = {rho_pt:.3f}')
print(f'Entropy vs published diff_pot: rho = {rho_dp:.3f}')

# ---------- 6. Update summary JSON ----------
summary_path = f'{RES}/palantir_summary.json'
with open(summary_path) as f:
    summary = json.load(f)

summary['validation_vs_published'] = {
    'pseudotime_vs_published': {
        'spearman_rho': float(rho_pt), 'p': float(p_pt), 'n': len(common),
    },
    'entropy_vs_diff_potential': {
        'spearman_rho': float(rho_dp), 'p': float(p_dp), 'n': len(common),
    },
    'branch_probabilities_vs_published': branch_rhos,
    'published_terminal_states': orig_ct,
    'mapping_note': 'Our "DC" maps to published "pDC + cDC" (sum). Our "Lymph" maps to published "CLP". Mega and Mk are not modeled in our 4-branch analysis.',
}

with open(summary_path, 'w') as f:
    json.dump(summary, f, indent=2)
print(f'\nUpdated {summary_path}')

# ---------- 7. Updated reproducibility figure ----------
n_branches = len(branch_rhos)
fig, axes = plt.subplots(2, max(3, n_branches), figsize=(4 * max(3, n_branches), 8))

# Row 0: pseudotime + entropy (same as before)
ax = axes[0, 0]
ax.scatter(orig_pt_aligned, our_aligned['palantir_pseudotime'],
           s=4, alpha=0.4, c='steelblue', linewidth=0)
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
ax.set_xlabel('Published pseudotime (Setty 2019)')
ax.set_ylabel('Re-derived pseudotime')
ax.set_title(f'Pseudotime (rho = {rho_pt:.3f})', fontweight='bold')

ax = axes[0, 1]
ax.scatter(orig_dp_aligned, our_aligned['palantir_entropy'],
           s=4, alpha=0.4, c='darkorange', linewidth=0)
ax.set_xlabel('Published diff. potential')
ax.set_ylabel('Re-derived entropy')
ax.set_title(f'Entropy (rho = {rho_dp:.3f})', fontweight='bold')

# Hide unused top-row axes
for j in range(2, axes.shape[1]):
    axes[0, j].axis('off')

# Row 1: per-branch
branch_colors = {'Ery': '#d62728', 'Mono': '#9467bd', 'DC': '#2ca02c', 'Lymph': '#1f77b4'}
for i, (label, stats) in enumerate(branch_rhos.items()):
    ax = axes[1, i]
    our_vals = our_aligned[f'fate_{label}'].values
    pub_vals = published_for_our_label[label]
    ax.scatter(pub_vals, our_vals, s=4, alpha=0.4,
               c=branch_colors.get(label, 'gray'), linewidth=0)
    ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
    map_to = stats['mapped_to_published']
    map_str = '+'.join(map_to) if isinstance(map_to, list) else map_to
    ax.set_xlabel(f'Published prob ({map_str})')
    ax.set_ylabel(f'Re-derived prob ({label})')
    ax.set_title(f'Branch {label} (rho = {stats["spearman_rho"]:.3f})', fontweight='bold')
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)

# Hide unused bottom-row axes
for j in range(n_branches, axes.shape[1]):
    axes[1, j].axis('off')

plt.suptitle('Reproducibility vs Setty 2019 published Palantir output',
             fontsize=13, fontweight='bold', y=1.00)
plt.tight_layout()
out_png = f'{FIG}/reproducibility_vs_published.png'
plt.savefig(out_png, dpi=300, bbox_inches='tight')
print(f'Saved figure: {out_png}')
plt.close()

print('\nDONE.')
