# Project 02 — Breast Cancer Tumor Microenvironment

**Tool:** Scanpy
**Dataset:** Wu et al. 2021 — single-cell atlas of human breast cancer
**Accession:** GEO **GSE176078** (also on CELLxGENE)
**Reference:** Wu SZ, Al-Eryani G, Roden DL, et al. *A single-cell and spatially resolved atlas of human breast cancers.* Nat Genet 53, 1334–1347 (2021). DOI: 10.1038/s41588-021-00911-1

## Why this project

Direct relevance to your manuscript domain (MDA-MB-231 breast cancer scRNA-seq) but on a fully public, well-cited dataset. Clients searching "breast cancer scRNA-seq" land on this.

## Goals
1. Build a multi-patient atlas covering tumor, stromal, and immune compartments
2. Identify epithelial subtypes (basal/luminal) and immune infiltrate composition
3. Visualize TME heterogeneity across patients with subtype information (HER2+, ER+, TNBC)

## Data download

```bash
# Option A: from CELLxGENE (preprocessed h5ad, recommended)
# Browse: https://cellxgene.cziscience.com/  → search "Wu 2021 breast"

# Option B: raw counts from GEO
wget -O data/GSE176078_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE176078&format=file"
```

## Estimated effort
**10–15 hours.** Multi-sample, larger dataset.
