#!/usr/bin/env python3
"""Quick probe: inspect Wu 2021 h5ad before writing the integration benchmark."""
import scanpy as sc

ad = sc.read('/home/marko-b2/upwork_portfolio/02_breast_cancer_TME/results/wu2021_3patient_TME_processed.h5ad')
print('Shape:', ad.shape)
print()
print('obs columns:')
for c in ad.obs.columns:
    n_unique = ad.obs[c].nunique()
    if n_unique < 30:
        print(f'  {c:<30} ({n_unique} unique): {sorted(ad.obs[c].astype(str).unique())[:10]}')
    else:
        print(f'  {c:<30} ({n_unique} unique, too many to list)')
print()
print('obsm keys:', list(ad.obsm.keys()))
print('layers keys:', list(ad.layers.keys()) if hasattr(ad, 'layers') else 'no layers')
print('uns keys:', list(ad.uns.keys())[:15])
print()

# X stats
import scipy.sparse as sp
X = ad.X
max_val = X.tocoo().data.max() if sp.issparse(X) else X.max()
print(f'X max value: {max_val:.3f}  (log-norm ~ 5-10, raw counts >> 100)')

# Check for raw counts somewhere
if ad.raw is not None:
    rX = ad.raw.X
    rmax = rX.tocoo().data.max() if sp.issparse(rX) else rX.max()
    print(f'raw.X max: {rmax:.1f}, shape: {ad.raw.shape}')
else:
    print('No .raw')

if hasattr(ad, 'layers') and 'counts' in ad.layers:
    cX = ad.layers['counts']
    cmax = cX.tocoo().data.max() if sp.issparse(cX) else cX.max()
    print(f'layers["counts"] max: {cmax:.1f}')
