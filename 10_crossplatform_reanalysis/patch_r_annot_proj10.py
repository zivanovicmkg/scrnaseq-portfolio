#!/usr/bin/env python3
"""Patch the R notebook to fix the annotation step that crashes on duplicate labels."""
import json

nb_path = '/home/marko-b2/upwork_portfolio/10_crossplatform_reanalysis/notebooks/02_r_analysis.ipynb'

with open(nb_path) as f:
    nb = json.load(f)

new_annot_src = """marker_signatures <- list(
  `CD4 T`        = c('IL7R','CCR7','TCF7','LTB'),
  `CD8 T`        = c('CD8A','CD8B','GZMK'),
  `NK`           = c('GNLY','NKG7','FCGR3A','GZMB','PRF1'),
  `B`            = c('MS4A1','CD79A','CD79B','IGHM','CD19'),
  `Mono CD14`    = c('CD14','LYZ','S100A8','S100A9'),
  `Mono FCGR3A`  = c('FCGR3A','MS4A7','CDKN1C'),
  `DC`           = c('FCER1A','CST3','CLEC10A'),
  `Platelet`     = c('PPBP','PF4','ITGA2B')
)

# Build cluster-to-label by max signature overlap
cluster_ids <- sort(as.integer(as.character(unique(markers_top$cluster))))
cluster_to_label <- character(length(cluster_ids))
names(cluster_to_label) <- as.character(cluster_ids)

for (cl in cluster_ids) {
  cl_markers <- markers_top %>% filter(cluster == cl) %>% head(30) %>% pull(gene)
  scores <- sapply(marker_signatures, function(sig) length(intersect(sig, cl_markers)))
  best <- names(which.max(scores))
  if (scores[best] == 0) best <- paste0('Unknown_', cl)
  cluster_to_label[as.character(cl)] <- best
}

# Disambiguate duplicates by appending cluster id
label_counts <- table(cluster_to_label)
dup_labels <- names(label_counts[label_counts > 1])
for (lab in dup_labels) {
  clusters_w_lab <- names(cluster_to_label)[cluster_to_label == lab]
  for (c in clusters_w_lab) {
    cluster_to_label[c] <- paste0(lab, '_', c)
  }
}

# Now map to obj$celltype via simple lookup (avoids RenameIdents pitfall)
obj$celltype <- factor(unname(cluster_to_label[as.character(obj$seurat_clusters)]))

cat('Final cluster -> celltype mapping:\\n')
for (cl in cluster_ids) {
  cat(sprintf('  %s: %s\\n', cl, cluster_to_label[as.character(cl)]))
}
cat('\\nCell-type counts:\\n')
print(table(obj$celltype))
"""

# Find and replace the annotation cell
patched = False
for c in nb['cells']:
    if c['cell_type'] == 'code':
        src = ''.join(c['source']) if isinstance(c['source'], list) else c['source']
        if 'marker_signatures' in src and 'cluster_to_label' in src and 'obj$celltype' in src:
            c['source'] = new_annot_src.splitlines(keepends=True)
            c['outputs'] = []
            c['execution_count'] = None
            patched = True
            print('Patched annotation cell')
            break

# Also patch the save cell to use simple list construction (older was OK but fragile)
new_save_src = """umap_coords <- Embeddings(obj, 'umap')
pca_coords <- Embeddings(obj, 'pca')

per_cell <- data.frame(
  cell = colnames(obj),
  seurat_cluster = as.character(obj$seurat_clusters),
  celltype = as.character(obj$celltype),
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2]
)
for (i in 1:min(30, ncol(pca_coords))) {
  per_cell[[paste0('PC', i)]] <- pca_coords[, i]
}
write.csv(per_cell, file.path(RES_DIR, 'r_per_cell.csv'), row.names = FALSE)

saveRDS(obj, file.path(RES_DIR, 'r_pbmc3k_processed.rds'))

summary_list <- list(
  ecosystem = 'Seurat (R)',
  seurat_version = as.character(packageVersion('Seurat')),
  n_cells_after_qc = ncol(obj),
  n_genes_after_qc = nrow(obj),
  n_hvg = length(VariableFeatures(obj)),
  n_pcs = 50,
  n_seurat_clusters = length(levels(obj$seurat_clusters)),
  cluster_sizes = as.list(table(obj$seurat_clusters)),
  cluster_to_celltype = as.list(cluster_to_label),
  celltypes = sort(unique(as.character(obj$celltype))),
  runtime_sec = list(
    qc = qc_time,
    preprocessing = prep_time,
    umap_clustering = umap_time,
    marker_genes = marker_time,
    total = qc_time + prep_time + umap_time + marker_time
  ),
  seed = 42
)
write_json(summary_list, file.path(RES_DIR, 'r_summary.json'), pretty = TRUE, auto_unbox = TRUE)

cat('Saved:\\n')
cat(' ', file.path(RES_DIR, 'r_per_cell.csv'), '\\n')
cat(' ', file.path(RES_DIR, 'r_pbmc3k_processed.rds'), '\\n')
cat(' ', file.path(RES_DIR, 'r_summary.json'), '\\n')
cat(' ', file.path(RES_DIR, 'r_markers_top50.csv'), '\\n')
cat(sprintf('\\nTotal pipeline time: %.1fs\\n', summary_list$runtime_sec$total))
"""

for c in nb['cells']:
    if c['cell_type'] == 'code':
        src = ''.join(c['source']) if isinstance(c['source'], list) else c['source']
        if 'r_summary.json' in src and 'r_per_cell.csv' in src:
            c['source'] = new_save_src.splitlines(keepends=True)
            c['outputs'] = []
            c['execution_count'] = None
            print('Patched save cell')
            break

# Clear all execution outputs to start clean
for c in nb['cells']:
    if c['cell_type'] == 'code':
        c['outputs'] = []
        c['execution_count'] = None

with open(nb_path, 'w') as f:
    json.dump(nb, f, indent=1)

print(f'Patched: {patched}')
