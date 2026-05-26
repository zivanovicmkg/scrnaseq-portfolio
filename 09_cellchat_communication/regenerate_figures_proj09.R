#!/usr/bin/env Rscript
# Project 09 - downstream figures from cached CellChat objects
# Loads cellchat_ERpos.rds, cellchat_HER2pos.rds, cellchat_TNBC.rds
# Generates everything from Section 5 onward (skips ~15 min of computeCommunProb)

suppressPackageStartupMessages({
  library(CellChat)
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(NMF)
  library(circlize)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)
set.seed(42)

BASE <- '/home/marko-b2/upwork_portfolio/09_cellchat_communication'
FIG  <- file.path(BASE, 'figures')
RES  <- file.path(BASE, 'results')
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, path, w = 8, h = 6, dpi = 300) {
  ggsave(path, plot = p, width = w, height = h, dpi = dpi, bg = 'white')
}

# ---------- 1. Load cached CellChat objects ----------
cat('Loading cached CellChat objects...\n')
subtypes <- c('ER+', 'HER2+', 'TNBC')
cellchat_list <- list(
  `ER+`   = readRDS(file.path(RES, 'cellchat_ERpos.rds')),
  `HER2+` = readRDS(file.path(RES, 'cellchat_HER2pos.rds')),
  `TNBC`  = readRDS(file.path(RES, 'cellchat_TNBC.rds'))
)
for (s in subtypes) {
  cc <- cellchat_list[[s]]
  cat(sprintf('  %-7s %d cells, %d cell types, %d L-R pairs, %d pathways\n',
              s, ncol(cc@data), length(levels(cc@idents)),
              nrow(cc@LR$LRsig), length(cc@netP$pathways)))
}

# ---------- 2. Hero circle plots (in case it failed) ----------
cat('\nSection 5a: hero circle plots...\n')
png(file.path(FIG, 'hero_circle_plots.png'), width = 16, height = 6, units = 'in', res = 300, bg = 'white')
par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))
for (s in subtypes) {
  cc <- cellchat_list[[s]]
  n_lr <- nrow(cc@LR$LRsig)
  n_path <- length(cc@netP$pathways)
  netVisual_circle(
    cc@net$weight,
    vertex.weight = as.numeric(table(cc@idents)),
    weight.scale = TRUE, label.edge = FALSE,
    title.name = sprintf('%s\n(%d L-R pairs, %d pathways)', s, n_lr, n_path)
  )
}
dev.off()
cat('  Saved hero_circle_plots.png\n')

# ---------- 3. Per-subtype heatmaps (separate PNGs, no merge) ----------
cat('\nSection 5b: per-subtype interaction heatmaps...\n')
for (s in subtypes) {
  cc <- cellchat_list[[s]]
  ht <- netVisual_heatmap(cc, measure = 'weight', color.heatmap = 'Reds', title.name = s)
  out_path <- file.path(FIG, paste0('interaction_heatmap_', gsub('[+]', 'pos', s), '.png'))
  png(out_path, width = 7, height = 6, units = 'in', res = 300, bg = 'white')
  ComplexHeatmap::draw(ht)
  dev.off()
  cat('  Saved', basename(out_path), '\n')
}

# ---------- 4. Lift to common labels + merge ----------
cat('\nSection 6: lift + merge for comparison...\n')
all_labels <- Reduce(union, lapply(cellchat_list, function(x) levels(x@idents)))
cat('  Union of all labels:', paste(all_labels, collapse = ', '), '\n')

for (s in subtypes) {
  cc <- cellchat_list[[s]]
  if (!all(all_labels %in% levels(cc@idents))) {
    cat('  Lifting', s, 'to common labels\n')
    cellchat_list[[s]] <- liftCellChat(cc, all_labels)
  }
}

merged <- mergeCellChat(cellchat_list, add.names = subtypes, cell.prefix = TRUE)
cat('  Merged object built.\n')
saveRDS(merged, file.path(RES, 'cellchat_merged.rds'))

# ---------- 5. Total interactions comparison ----------
cat('\nSection 7: total interactions per subtype...\n')
gg1 <- compareInteractions(merged, show.legend = FALSE, group = seq_along(subtypes),
                           color.use = c('#1f77b4','#ff7f0e','#d62728'))
gg2 <- compareInteractions(merged, show.legend = FALSE, group = seq_along(subtypes),
                           measure = 'weight', color.use = c('#1f77b4','#ff7f0e','#d62728'))
p <- gg1 + gg2 + plot_annotation(title = 'Global signaling activity per subtype',
                                  theme = theme(plot.title = element_text(face = 'bold', size = 13)))
save_fig(p, file.path(FIG, 'total_interactions_per_subtype.png'), w = 9, h = 4.5)
cat('  Saved total_interactions_per_subtype.png\n')

# ---------- 6. Differential interaction circles (pairwise) ----------
cat('\nSection 8: pairwise differential interaction circles...\n')
pairwise <- list(c('ER+','TNBC'), c('ER+','HER2+'), c('HER2+','TNBC'))
png(file.path(FIG, 'differential_interaction_circles.png'),
    width = 16, height = 6, units = 'in', res = 300, bg = 'white')
par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))
for (pair in pairwise) {
  m_pair <- mergeCellChat(cellchat_list[pair], add.names = pair, cell.prefix = TRUE)
  netVisual_diffInteraction(m_pair, weight.scale = TRUE, measure = 'weight',
                            title.name = sprintf('%s vs %s (red: up in %s)', pair[2], pair[1], pair[2]))
}
dev.off()
cat('  Saved differential_interaction_circles.png\n')

# ---------- 7. Pathway-level rankNet ----------
cat('\nSection 9: pathway rankings...\n')
rank_stacked <- rankNet(merged, mode = 'comparison', stacked = TRUE, do.stat = TRUE,
                        color.use = c('#1f77b4','#ff7f0e','#d62728'))
save_fig(rank_stacked, file.path(FIG, 'pathway_rank_stacked.png'), w = 8, h = 12)
cat('  Saved pathway_rank_stacked.png\n')

rank_unstacked <- rankNet(merged, mode = 'comparison', stacked = FALSE, do.stat = TRUE,
                          color.use = c('#1f77b4','#ff7f0e','#d62728'))
save_fig(rank_unstacked, file.path(FIG, 'pathway_rank_unstacked.png'), w = 12, h = 7)
cat('  Saved pathway_rank_unstacked.png\n')

# ---------- 8. Top L-R pairs bubble plot ----------
cat('\nSection 10: top L-R pairs bubble plot + per-subtype CSVs...\n')
df_all <- list()
for (s in subtypes) {
  cc <- cellchat_list[[s]]
  df_s <- subsetCommunication(cc)
  if (nrow(df_s) > 0) {
    df_s <- df_s[order(-df_s$prob), ]
    write.csv(df_s, file.path(RES, paste0('LR_pairs_', gsub('[+]', 'pos', s), '.csv')), row.names = FALSE)
    df_s$subtype <- s
    df_all[[s]] <- df_s
  }
}
df_combined <- do.call(rbind, df_all)

top_lr <- df_combined %>%
  group_by(interaction_name) %>%
  summarise(total_prob = sum(prob, na.rm = TRUE), n_subtypes = n_distinct(subtype)) %>%
  arrange(desc(total_prob)) %>% head(25) %>% pull(interaction_name)

df_top <- df_combined %>% filter(interaction_name %in% top_lr) %>%
  mutate(pair = paste(source, '->', target),
         interaction_name = factor(interaction_name, levels = rev(top_lr)))

p_bubble <- ggplot(df_top, aes(x = subtype, y = interaction_name,
                                size = prob, color = -log10(pval + 1e-3))) +
  geom_point() +
  scale_color_gradient(low = 'steelblue', high = 'firebrick', name = '-log10(p)') +
  scale_size_continuous(range = c(1, 6), name = 'Comm prob') +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(color = 'gray90', size = 0.3)) +
  labs(title = 'Top 25 L-R interactions across subtypes', x = NULL, y = NULL)
save_fig(p_bubble, file.path(FIG, 'top_LR_bubble.png'), w = 9, h = 10)
cat('  Saved top_LR_bubble.png\n')

# ---------- 9. Outgoing signaling heatmaps (per subtype, separate PNGs) ----------
cat('\nSection 11: outgoing signaling heatmaps (separate PNGs per subtype)...\n')
for (s in subtypes) {
  cc <- cellchat_list[[s]]
  ht <- netAnalysis_signalingRole_heatmap(
    cc, pattern = 'outgoing', height = 8, width = 5,
    color.heatmap = 'OrRd', title = paste(s, '- outgoing'))
  out_path <- file.path(FIG, paste0('signaling_role_outgoing_', gsub('[+]', 'pos', s), '.png'))
  png(out_path, width = 6, height = 7, units = 'in', res = 300, bg = 'white')
  ComplexHeatmap::draw(ht)
  dev.off()
  cat('  Saved', basename(out_path), '\n')
}

# ---------- 10. Summary JSON ----------
cat('\nSection 12: summary JSON...\n')
summary_list <- list(
  dataset = 'Wu 2021 breast cancer, 3 patients, 12,962 cells',
  subtypes = subtypes,
  cellchat_version = as.character(packageVersion('CellChat')),
  per_subtype = list()
)

for (s in subtypes) {
  cc <- cellchat_list[[s]]
  summary_list$per_subtype[[s]] <- list(
    n_cells = ncol(cc@data),
    n_celltypes = length(levels(cc@idents)),
    n_significant_LR_pairs = nrow(cc@LR$LRsig),
    n_signaling_pathways = length(cc@netP$pathways),
    total_n_interactions = sum(cc@net$count),
    total_interaction_strength = round(sum(cc@net$weight), 4),
    top_pathways = head(cc@netP$pathways, 10)
  )
}

write_json(summary_list, file.path(RES, 'cellchat_summary.json'),
           pretty = TRUE, auto_unbox = TRUE)

cat('\n=== SUMMARY ===\n')
for (s in subtypes) {
  ps <- summary_list$per_subtype[[s]]
  cat(sprintf('\n--- %s ---\n', s))
  cat(sprintf('  Cells:              %d\n', ps$n_cells))
  cat(sprintf('  Cell types:         %d\n', ps$n_celltypes))
  cat(sprintf('  Significant L-R:    %d\n', ps$n_significant_LR_pairs))
  cat(sprintf('  Pathways:           %d\n', ps$n_signaling_pathways))
  cat(sprintf('  Total interactions: %d\n', ps$total_n_interactions))
  cat(sprintf('  Total strength:     %.4f\n', ps$total_interaction_strength))
  cat(sprintf('  Top 10 pathways:    %s\n', paste(ps$top_pathways, collapse = ', ')))
}

cat('\nDONE.\n')
