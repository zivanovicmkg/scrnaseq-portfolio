# Project 08 — Multi-Sample Integration & Batch Correction

**Tool:** Scanpy + Harmony + BBKNN + Scanorama (benchmark)
**Dataset:** Re-uses multi-patient data from Project 02 (Wu 2021 breast cancer) or Project 05 (CRC atlas)

## Why this project

Pure **technical showcase**. Batch effects are the #1 problem clients ask about when sending multi-sample data. This entry says: "I know how to integrate, and I know which method to choose."

## Goals
1. Demonstrate the integration problem (uncorrected UMAP showing batch separation)
2. Apply 3 integration methods: Harmony, BBKNN, Scanorama
3. Quantitative benchmark: kBET, LISI, silhouette by batch vs. cell type
4. Recommend a method based on observed metrics

## Estimated effort
**6–8 hours** (reuses preprocessed data from earlier projects).
