# Upwork portfolio briefs — ready-to-paste text blocks

Each block is under 600 characters and matches the structure that worked for the niraparib entry: 1-sentence framing + methods + deliverables.

---

## 01 — PBMC 3K: Reference scRNA-seq Pipeline

> End-to-end single-cell RNA-seq workflow on the 10X Genomics PBMC 3K dataset, demonstrating a clean, reproducible pipeline on a benchmark immune-cell reference.
>
> Full Python (Scanpy) pipeline: QC, normalization, highly variable gene selection, PCA, UMAP embedding, Leiden clustering, marker-based annotation of T cells, B cells, NK, monocytes, and dendritic cells. Deliverables: annotated UMAP, marker heatmap, top-gene dotplots, and a reproducible notebook (300 dpi figures, journal-ready).

---

## 02 — Breast Cancer Tumor Microenvironment

> Single-cell atlas analysis of human breast cancer (Wu et al. 2021, GSE176078) covering tumor, stromal, and immune compartments across multiple patients.
>
> Scanpy pipeline: per-sample QC, doublet removal, multi-sample integration, Leiden clustering at multiple resolutions, marker-based annotation, and compositional comparison across breast cancer subtypes (HER2+, ER+, TNBC). Deliverables: integrated UMAP, subtype-stratified composition plots, marker heatmaps, and a publication-style mini-report.

---

## 03 — Glioblastoma Cellular States

> Reproduction of the four canonical glioblastoma cellular states (Mesenchymal, Astrocyte-like, OPC-like, NPC-like) on Neftel et al. 2019 (GSE131928) using Seurat.
>
> Full R/Seurat workflow: QC, SCTransform normalization, integration, clustering, and module-score signatures for each state. Plotted the canonical 2D state diagram and characterized state plasticity per tumor. Deliverables: state-scored UMAPs, butterfly state plot, marker heatmaps, and a reproducible Rmarkdown report.

---

## 04 — Melanoma Immune Landscape

> Re-analysis of Jerby-Arnon et al. 2018 (GSE115978) characterizing malignant and immune compartments in melanoma, with a focus on T cell exclusion and checkpoint-resistance signatures.
>
> Seurat pipeline: malignant vs non-malignant classification, T cell exhaustion scoring, and reconstruction of the published "T cell exclusion" program. Deliverables: annotated UMAPs, exhaustion-marker dotplots, exclusion-program score visualization, and a written interpretation suitable for immuno-oncology project briefs.

---

## 05 — Colorectal Cancer Single-Cell Atlas

> Multi-patient single-cell analysis of colorectal cancer (Pelka et al. 2021 / Lee et al. 2020), integrating epithelial, stromal, and immune compartments across dozens of tumors.
>
> Scanpy pipeline: scalable per-sample QC, doublet removal, Harmony integration across donors, Leiden clustering, and lineage annotation. Compositional analysis stratified by available clinical metadata. Deliverables: integrated atlas UMAP, lineage composition plots, marker dotplots, and a written summary of TME structure.

---

## 06 — RNA Velocity: Pancreatic Endocrinogenesis

> RNA velocity analysis of mouse pancreatic endocrinogenesis (Bastidas-Ponce et al. 2019) — the reference dataset for the scVelo dynamical model.
>
> Python (scVelo) pipeline: steady-state and dynamical velocity computation, latent time estimation, velocity streams overlaid on UMAP, driver-gene identification, and PAGA-velocity graph for endocrine lineage trajectories (alpha, beta, delta, epsilon). Deliverables: velocity-stream UMAP, latent-time map, driver-gene plots, and a methodological report on velocity interpretation.

---

## 07 — Trajectory Inference: Hematopoiesis

> Trajectory and fate-probability analysis of human CD34+ HSPCs (Setty et al. 2019) using Palantir, complementing velocity-based approaches with a probabilistic trajectory framework.
>
> Python pipeline (Palantir + Scanpy): diffusion-map embedding, pseudotime computation, fate-probability estimation toward erythroid, lymphoid, and myeloid lineages, and gene-trend analysis along each branch. Deliverables: pseudotime UMAP, fate-probability heatmap, lineage-specific gene trends, and a comparative note on trajectory vs velocity methods.

---

## 08 — Multi-Sample Integration & Batch Correction

> Benchmark of three single-cell integration methods (Harmony, BBKNN, Scanorama) on a multi-patient cancer dataset, with quantitative evaluation of batch removal versus biological signal preservation.
>
> Scanpy pipeline: uncorrected baseline, integration with all three methods, and quantitative scoring using kBET, LISI, and silhouette metrics (batch vs. cell type). Deliverables: side-by-side UMAP comparison, integration metric tables, and a recommendation report on when to choose which method.

---

## 09 — Cell-Cell Communication Analysis

> Ligand-receptor signaling analysis of a breast cancer tumor microenvironment, inferring intercellular communication across tumor, stromal, and immune cell types using CellChat.
>
> R/CellChat pipeline: signaling network inference, dominant ligand-receptor pair identification, pathway-level aggregation, and visualization with chord, circle, and heatmap plots. Deliverables: cell-cell communication network plots, top signaling pathway summaries, and a written interpretation of dominant TME interactions.

---

## 10 — Cross-Platform Re-Analysis & Benchmarking

> Parallel re-analysis of a published scRNA-seq dataset in both Scanpy (Python) and Seurat (R), comparing clustering, annotation, and biological conclusions across two independent ecosystems.
>
> Workflow: identical preprocessing standards applied in both tools, side-by-side clustering, marker gene overlap analysis, and reconciliation of differences. Deliverables: cluster-membership comparison, marker overlap tables, dual UMAPs, and a transparent report on where the two pipelines agree and where methodological choices change conclusions.
