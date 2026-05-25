# Project 06 — RNA Velocity: Pancreatic Endocrinogenesis

**Tool:** scVelo (Python)
**Dataset:** Bastidas-Ponce et al. 2019 — mouse pancreatic endocrinogenesis (E15.5)
**Reference:** Bastidas-Ponce A, Tritschler S, Dony L, et al. *Comprehensive single cell mRNA profiling reveals a detailed roadmap for pancreatic endocrinogenesis.* Development 146 (2019). DOI: 10.1242/dev.173849

## Why this project

Dedicated **methodological showcase** for RNA velocity. The scVelo reference dataset — every velocity analysis you see in papers compares to this. Demonstrates command of dynamical model, latent time, and PAGA-velocity integration.

## Goals
1. Compute steady-state and dynamical velocity models
2. Plot velocity streams overlaid on UMAP (the "hero" RNA velocity figure)
3. Calculate latent time and rank driver genes
4. Project trajectories from progenitor to differentiated endocrine cells (alpha, beta, delta, epsilon)

## Data download
```python
import scvelo as scv
adata = scv.datasets.pancreas()
```
Built-in to scVelo — instant access.

## Estimated effort
**8–10 hours.** Smaller dataset, but velocity computation is the visual centerpiece.
