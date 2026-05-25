# Project 01 — PBMC 3K: Reference scRNA-seq Pipeline

**Tool:** Scanpy
**Dataset:** 10X Genomics PBMC 3K (3000 peripheral blood mononuclear cells from a healthy donor)
**Reference:** 10X Genomics public dataset; benchmark used by Wolf et al. 2018 (Genome Biology)

## Why this project

A clean, end-to-end reference workflow. Demonstrates command of the standard scRNA-seq pipeline on a well-known dataset every bioinformatician recognizes. Acts as the "competence proof" entry in the portfolio.

## Goals
1. Demonstrate full QC → clustering → annotation pipeline
2. Produce publication-grade UMAP, violin, and dotplot figures
3. Annotate canonical immune cell types (T cells, B cells, NK, monocytes, dendritic cells)

## Data download

Inside the notebook:
```python
adata = sc.datasets.pbmc3k()
```
No manual download needed — bundled with Scanpy.

## Estimated effort
**4–6 hours.** First project; also a sanity check that the environment works.
