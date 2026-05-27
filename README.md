# scRNA-Seq Portfolio — Marko Živanović

End-to-end single-cell RNA-seq analyses on public datasets, prepared as a portfolio of reproducible workflows across cancer biology, immunology, development, and methodological showcases. Every project ships with: figures, a Jupyter or R notebook, and a short PDF report.

**Tools:** Scanpy (Python), Seurat (R), scVelo, Palantir, Harmony, BBKNN, Scanorama, scVI, CellChat, scib-metrics
**Reproducibility:** conda environments, pinned versions, public data only
**Author:** Marko Živanović — PhD, Principal Research Fellow, BioIRC / University of Kragujevac

---

## Project index

| #  | Title | Tool | Dataset | Status |
|----|-------|------|---------|--------|
| 01 | PBMC 3K — Reference scRNA-seq Pipeline | Scanpy | 10X PBMC 3K | ✅ |
| 02 | Breast Cancer Tumor Microenvironment | Scanpy | Wu et al. 2021, GSE176078 | ✅ |
| 03 | Glioblastoma Cellular States | Seurat | Neftel et al. 2019, GSE131928 | ✅ |
| 04 | Melanoma Immune Landscape | Seurat | Jerby-Arnon et al. 2018, GSE115978 | ✅ |
| 05 | Colorectal Cancer Single-Cell Atlas | Scanpy | Pelka et al. 2021 / Lee et al. 2020 | ✅ |
| 06 | RNA Velocity — Pancreatic Endocrinogenesis | scVelo | Bastidas-Ponce et al. 2019 | ✅ |
| 07 | Trajectory Inference — Hematopoiesis | Palantir + Scanpy | Setty et al. 2019 | ✅ |
| 08 | Multi-Sample Integration & Batch Correction | Scanpy + Harmony | Re-using #02 or #05 | ✅ |
| 09 | Cell-Cell Communication Analysis | CellChat (R) | Re-using #02 (TME) | ✅ |
| 10 | Cross-Platform Re-Analysis & Benchmarking | Scanpy + Seurat | TBD (published set) | ✅ |

Tick the box (`☐` → `☑`) when each project is fully shipped (figures + notebook + report).

---

## Repository layout

```
upwork_portfolio/
├── README.md                    # this file
├── environment.yml              # conda env (Python + R)
├── requirements.txt             # pip-only alternative
├── Dockerfile                   # containerized alternative
├── .gitignore
├── _shared/
│   ├── scripts/                 # helper scripts reused across projects
│   ├── styles/                  # matplotlib/seaborn style for consistent figures
│   └── templates/               # notebook + report skeletons
├── docs/
│   └── portfolio_brief.md       # short text descriptions ready for Upwork
└── NN_project_name/
    ├── README.md                # project-specific brief
    ├── data/                    # raw + processed data (gitignored, large)
    ├── notebooks/               # Jupyter / Rmd analysis notebooks
    ├── figures/                 # PNG + TIFF outputs (300 dpi)
    ├── results/                 # exported tables, h5ad/rds objects
    └── report/                  # PDF mini-report
```

---

## How to run

```bash
# create the main environment
mamba env create -f environment.yml
mamba activate scportfolio

# verify
python -c "import scanpy; print(scanpy.__version__)"
Rscript -e 'library(Seurat); packageVersion("Seurat")'
```

---

## Portfolio-ready output checklist

Each project must produce, before being marked complete:

- [ ] **Hero figure** — 1 visually striking publication-ready figure (UMAP, heatmap, velocity, etc.) saved as PNG (1200×900) for Upwork thumbnail + TIFF 300 dpi for manuscript-grade
- [ ] **Supporting figures** — 3–6 additional analytical figures
- [ ] **Notebook** — clean, narrative-driven, runnable end-to-end
- [ ] **Mini-report** — 2–3 page PDF (introduction → methods → results → biological interpretation)
- [ ] **Project README** — short summary + how-to-run + figure preview
- [ ] **Public link** — committed and pushed to GitHub

---

## License & data

All datasets are public (GEO, Zenodo, Single Cell Portal, 10X Genomics) — proper citations are included in each project's README. Code in this repository is released under MIT.
