#!/usr/bin/env Rscript
# Additional analyses for Project 04:
#   B) Relax per-tumor cell threshold from 20 to 10 malignant cells
#   C) Add per-tumor exclusion distribution stratified by T cell infiltration

suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(ggplot2)
})

set.seed(42)

# Paths (run from notebooks folder)
DATA_DIR    <- '../data'
FIG_DIR     <- '../figures'
RESULTS_DIR <- '../results'

# Same save_fig helper used elsewhere
save_fig <- function(plot_obj, name, w = 6, h = 5, dpi = 300, out_dir = FIG_DIR) {
    for (fmt in c('png', 'tiff', 'pdf')) {
        path <- file.path(out_dir, paste0(name, '.', fmt))
        if (fmt == 'tiff') {
            ggsave(path, plot_obj, width = w, height = h, dpi = dpi,
                   device = 'tiff', compression = 'lzw')
        } else {
            ggsave(path, plot_obj, width = w, height = h, dpi = dpi,
                   device = fmt)
        }
        cat('  saved', path, '\n')
    }
}

# 1. Load processed object
cat('Loading processed Seurat object...\n')
obj <- readRDS(file.path(RESULTS_DIR, 'jerby_arnon_processed.rds'))
cat('Loaded:', ncol(obj), 'cells\n')

# 2. Detect sample column
sample_col <- if ('samples' %in% colnames(obj[[]])) 'samples' else 'orig.ident'

df_tumor <- obj[[]] %>%
    as_tibble(rownames = 'cell') %>%
    rename(sample = !!sample_col)

# ----- B) Per-tumor with RELAXED threshold (>= 10 malignant cells) -----
cat('\n=== B) Relaxed threshold per-tumor correlation ===\n')

per_sample_relaxed <- df_tumor %>%
    group_by(sample) %>%
    summarise(
        total_cells    = n(),
        n_malignant    = sum(lineage == 'Malignant'),
        n_tcell        = sum(lineage == 'T_cell'),
        tcell_fraction = n_tcell / total_cells,
        mean_exclusion = mean(ExclusionScore[lineage == 'Malignant'], na.rm = TRUE),
        .groups = 'drop'
    ) %>%
    filter(n_malignant >= 10, total_cells >= 30)

cat(sprintf('Tumors with >=10 malignant cells: %d\n', nrow(per_sample_relaxed)))

if (nrow(per_sample_relaxed) >= 4) {
    corr_relaxed <- cor.test(per_sample_relaxed$mean_exclusion,
                              per_sample_relaxed$tcell_fraction,
                              method = 'pearson')
    cat(sprintf('Pearson r = %.3f, p = %.4f\n',
                corr_relaxed$estimate, corr_relaxed$p.value))

    # Spearman too — robust to outliers, often more appropriate for small n
    corr_spearman <- cor.test(per_sample_relaxed$mean_exclusion,
                               per_sample_relaxed$tcell_fraction,
                               method = 'spearman', exact = FALSE)
    cat(sprintf('Spearman rho = %.3f, p = %.4f\n',
                corr_spearman$estimate, corr_spearman$p.value))

    p_corr <- ggplot(per_sample_relaxed, aes(x = mean_exclusion, y = tcell_fraction)) +
        geom_point(aes(size = total_cells), color = '#1f77b4', alpha = 0.7) +
        geom_smooth(method = 'lm', se = TRUE, color = '#d62728', fill = '#fcd0cd') +
        scale_size_continuous(name = 'Cells\nper tumor', range = c(2, 7)) +
        labs(
            title = 'T cell exclusion program vs T cell infiltration',
            subtitle = sprintf('Per-tumor: Pearson r = %.2f (p = %.3f); Spearman rho = %.2f (p = %.3f); n = %d tumors',
                                corr_relaxed$estimate, corr_relaxed$p.value,
                                corr_spearman$estimate, corr_spearman$p.value,
                                nrow(per_sample_relaxed)),
            x = 'Mean exclusion score (malignant cells)',
            y = 'T cell fraction (all cells)'
        ) +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = 'bold'),
              plot.subtitle = element_text(size = 9))
    save_fig(p_corr, '06_exclusion_vs_tcell_infiltration_relaxed', w = 9, h = 6)
    write.csv(per_sample_relaxed,
              file.path(RESULTS_DIR, 'per_tumor_exclusion_tcell_relaxed.csv'),
              row.names = FALSE)
}

# ----- C) Stratified analysis: malignant-cell exclusion distribution by tumor T cell infiltration -----
cat('\n=== C) Stratified analysis: malignant cells exclusion by tumor T-cell infiltration ===\n')

# Use median tcell_fraction across the relaxed set as split point
tcell_split <- median(per_sample_relaxed$tcell_fraction)
high_inf_tumors <- per_sample_relaxed$sample[per_sample_relaxed$tcell_fraction >= tcell_split]
low_inf_tumors  <- per_sample_relaxed$sample[per_sample_relaxed$tcell_fraction <  tcell_split]

cat(sprintf('Median T cell fraction: %.3f\n', tcell_split))
cat(sprintf('High-infiltration tumors (>= median): %d\n', length(high_inf_tumors)))
cat(sprintf('Low-infiltration tumors  (<  median): %d\n', length(low_inf_tumors)))

# Build per-cell dataframe of MALIGNANT cells only, with the tumor infiltration class
mal_cells <- df_tumor %>%
    filter(lineage == 'Malignant', sample %in% per_sample_relaxed$sample) %>%
    mutate(
        tumor_infiltration = case_when(
            sample %in% high_inf_tumors ~ 'High T-cell infiltration',
            sample %in% low_inf_tumors  ~ 'Low T-cell infiltration',
            TRUE ~ NA_character_
        ),
        tumor_infiltration = factor(tumor_infiltration,
                                     levels = c('Low T-cell infiltration',
                                                'High T-cell infiltration'))
    )

cat(sprintf('Malignant cells in analysis: %d\n', nrow(mal_cells)))

# Statistical test: are exclusion scores systematically higher in low-infiltration tumors?
wt <- wilcox.test(ExclusionScore ~ tumor_infiltration, data = mal_cells)
cat(sprintf('\nWilcoxon rank-sum test (Low vs High infiltration tumors):\n'))
cat(sprintf('  W = %.0f, p = %.4g\n', wt$statistic, wt$p.value))

# Per-group summary
group_summary <- mal_cells %>%
    group_by(tumor_infiltration) %>%
    summarise(
        n_cells = n(),
        median_excl = median(ExclusionScore),
        mean_excl = mean(ExclusionScore),
        .groups = 'drop'
    )
cat('\nGroup summary:\n')
print(group_summary)

# Plot: boxplot + violin
p_strat <- ggplot(mal_cells, aes(x = tumor_infiltration, y = ExclusionScore,
                                   fill = tumor_infiltration)) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = c('Low T-cell infiltration' = '#d62728',
                                  'High T-cell infiltration' = '#1f77b4'),
                       guide = 'none') +
    labs(
        title = 'Malignant-cell exclusion score by tumor T cell infiltration',
        subtitle = sprintf('Wilcoxon p = %.2e ; n_low = %d cells from %d tumors, n_high = %d cells from %d tumors',
                            wt$p.value,
                            sum(mal_cells$tumor_infiltration == 'Low T-cell infiltration'),
                            length(low_inf_tumors),
                            sum(mal_cells$tumor_infiltration == 'High T-cell infiltration'),
                            length(high_inf_tumors)),
        x = NULL, y = 'Exclusion score (per malignant cell)'
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = 'bold'),
          plot.subtitle = element_text(size = 9))

save_fig(p_strat, '07_exclusion_stratified_by_infiltration', w = 7, h = 5.5)
write.csv(group_summary,
          file.path(RESULTS_DIR, 'exclusion_stratified_summary.csv'),
          row.names = FALSE)

cat('\nDone.\n')
