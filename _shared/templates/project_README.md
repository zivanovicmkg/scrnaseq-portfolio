# Project NN — {{TITLE}}

**Tool:** {{Scanpy | Seurat | scVelo | ...}}
**Dataset:** {{Dataset name + accession (GEO / Zenodo / 10X)}}
**Reference:** {{Author et al., Journal Year, DOI}}

---

## Summary

One paragraph (3–5 sentences) describing the biological context, what the analysis aims to demonstrate, and the headline result.

## Methods

Bulleted list of the actual steps performed (QC thresholds, normalization, integration method, clustering resolution, etc.).

- Quality control: ...
- Normalization: ...
- Feature selection: ...
- Dimensionality reduction: ...
- Clustering: ...
- Annotation: ...
- Downstream analysis: ...

## Key figures

| File | Description |
|------|-------------|
| `figures/01_umap_clusters.png`  | UMAP colored by Leiden clusters |
| `figures/02_marker_heatmap.png` | Top-marker heatmap per cluster |
| `figures/03_xxx.png`            | ... |

## How to run

```bash
mamba activate scportfolio
cd notebooks
jupyter lab analysis.ipynb     # or: jupyter nbconvert --execute analysis.ipynb
```

Data is downloaded inside the notebook from public sources; nothing private.

## Files

- `notebooks/analysis.ipynb` — main analysis
- `report/report.pdf` — 2–3 page interpretation
- `figures/` — PNG + TIFF outputs
- `results/` — exported tables, processed h5ad/rds

## Citation

If you reuse this analysis, please cite the original dataset publication and this repository.
