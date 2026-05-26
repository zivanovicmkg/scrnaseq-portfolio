#!/usr/bin/env Rscript
# Patch: regenerate pathway_rank plots showing all 3 subtypes
# The previous rankNet output only showed 2 colors (ER+/HER2+).
# Use comparison='pairwise' for each pair, OR plot the raw pathway info-flow matrix directly.

suppressPackageStartupMessages({
  library(CellChat)
  library(tidyverse)
  library(patchwork)
})

set.seed(42)

BASE <- '/home/marko-b2/upwork_portfolio/09_cellchat_communication'
FIG  <- file.path(BASE, 'figures')
RES  <- file.path(BASE, 'results')

subtypes <- c('ER+', 'HER2+', 'TNBC')
cellchat_list <- list(
  `ER+`   = readRDS(file.path(RES, 'cellchat_ERpos.rds')),
  `HER2+` = readRDS(file.path(RES, 'cellchat_HER2pos.rds')),
  `TNBC`  = readRDS(file.path(RES, 'cellchat_TNBC.rds'))
)

# Lift all to common labels
all_labels <- Reduce(union, lapply(cellchat_list, function(x) levels(x@idents)))
for (s in subtypes) {
  cc <- cellchat_list[[s]]
  if (!all(all_labels %in% levels(cc@idents))) {
    cellchat_list[[s]] <- liftCellChat(cc, all_labels)
  }
}

# Build per-pathway info-flow matrix manually:
# For each subtype, sum communication probability across all cell-type pairs per pathway
pathway_flow <- list()
for (s in subtypes) {
  cc <- cellchat_list[[s]]
  prob <- cc@netP$prob  # [source, target, pathway]
  if (length(dim(prob)) == 3) {
    pw_sum <- apply(prob, 3, function(m) sum(m, na.rm = TRUE))
    pathway_flow[[s]] <- data.frame(pathway = dimnames(prob)[[3]], info_flow = pw_sum, subtype = s)
  }
}
df_flow <- do.call(rbind, pathway_flow)
rownames(df_flow) <- NULL

# Wide format for sorting and side-by-side
df_wide <- df_flow %>% pivot_wider(names_from = subtype, values_from = info_flow, values_fill = 0)
df_wide$total <- df_wide$`ER+` + df_wide$`HER2+` + df_wide$`TNBC`
df_wide <- df_wide %>% arrange(desc(total))

cat('Top 30 pathways by total info flow across subtypes:\n')
print(head(df_wide, 30))

# Plot 1: stacked horizontal bars, top 40 pathways
top_n <- 40
df_top <- df_flow %>% filter(pathway %in% df_wide$pathway[1:top_n]) %>%
  mutate(pathway = factor(pathway, levels = rev(df_wide$pathway[1:top_n])),
         subtype = factor(subtype, levels = subtypes))

# Normalize per pathway so each row sums to 1 (relative information flow)
df_top_norm <- df_top %>%
  group_by(pathway) %>%
  mutate(rel_flow = info_flow / sum(info_flow)) %>%
  ungroup()

p_stacked <- ggplot(df_top_norm, aes(x = rel_flow, y = pathway, fill = subtype)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = c('ER+' = '#1f77b4', 'HER2+' = '#ff7f0e', 'TNBC' = '#d62728')) +
  geom_vline(xintercept = 1/3, linetype = 'dashed', color = 'gray40', size = 0.4) +
  labs(x = 'Relative information flow (proportion by subtype)',
       y = NULL,
       fill = NULL,
       title = paste('Top', top_n, 'signaling pathways: relative contribution per subtype')) +
  theme_minimal(base_size = 10) +
  theme(legend.position = 'top',
        panel.grid.major.y = element_blank())

ggsave(file.path(FIG, 'pathway_rank_stacked.png'), p_stacked, width = 9, height = 12, dpi = 300, bg = 'white')
cat('Saved pathway_rank_stacked.png\n')

# Plot 2: absolute info flow (side-by-side bars), top 25
df_top25 <- df_flow %>% filter(pathway %in% df_wide$pathway[1:25]) %>%
  mutate(pathway = factor(pathway, levels = rev(df_wide$pathway[1:25])),
         subtype = factor(subtype, levels = subtypes))

p_abs <- ggplot(df_top25, aes(x = info_flow, y = pathway, fill = subtype)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c('ER+' = '#1f77b4', 'HER2+' = '#ff7f0e', 'TNBC' = '#d62728')) +
  labs(x = 'Information flow (sum of communication probability)',
       y = NULL, fill = NULL,
       title = 'Top 25 pathways by total signaling activity') +
  theme_minimal(base_size = 10) +
  theme(legend.position = 'top')

ggsave(file.path(FIG, 'pathway_rank_unstacked.png'), p_abs, width = 11, height = 8, dpi = 300, bg = 'white')
cat('Saved pathway_rank_unstacked.png\n')

# Save data
write.csv(df_wide, file.path(RES, 'pathway_info_flow.csv'), row.names = FALSE)
cat('Saved pathway_info_flow.csv\n')

cat('\nDONE.\n')
