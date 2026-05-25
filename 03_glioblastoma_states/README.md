# Project 03 — Glioblastoma Cellular States

**Tool:** Seurat (R)
**Dataset:** Neftel et al. 2019 — IDH-wild-type glioblastoma
**Accession:** GEO **GSE131928** (Smart-seq2) + 10X portion
**Reference:** Neftel C, Laffy J, Filbin MG, et al. *An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma.* Cell 178, 835–849 (2019). DOI: 10.1016/j.cell.2019.06.024

## Why this project

The Seurat showcase. Reproduces the four canonical GBM cellular states (MES, AC, OPC, NPC) and demonstrates state-scoring methodology. Clients who work in R/Seurat will look at this entry first.

## Goals
1. Reproduce the four-state classification (Mesenchymal, Astrocyte-like, OPC-like, NPC-like)
2. Compute state scores and plot the canonical 2D state representation (butterfly plot)
3. Demonstrate Seurat's AddModuleScore and CellCycleScoring workflows

## Data download
```r
# From GEO
download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131928/suppl/GSE131928_RAW.tar",
  "data/GSE131928_RAW.tar"
)
```

## Estimated effort
**12–18 hours.** Seurat workflow + state scoring + signature 2D plot.
