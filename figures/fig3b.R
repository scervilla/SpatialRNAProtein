library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)

path_output <- "figures/figure3/"

datasets <- c("LUAD", "mCRPC")
metric_labels <- c(detection_kappa = "Cohen's Kappa",
                   spearman_all = "Spearman cor",
                   spearman_pos = "Spearman cor (pos)",
                   lee_L = "Lee's L")
metric_colours <- list(`Cohen's Kappa` = c(low = "#9ECAE1", high = "#08519C"),
                       `Spearman cor` = c(low = "#FCBBA1", high = "#A50F15"),
                       `Spearman cor (pos)` = c(low = "#A1D99B", high = "#006D2C"),
                       `Lee's L` = c(low = "#BCBDDC", high = "#54278F"))

# Per-sample concordance metrics (from concordance_metrics.R)
results <- bind_rows(lapply(datasets, function(dataset) {
  samples <- readLines(paste0("utils/samples_", dataset, ".txt"))
  bind_rows(lapply(samples, function(sample) {
    read.table(paste0("output/concordance/", dataset, "/", sample, "_summary.txt"))
  })) %>%
    mutate(dataset = dataset)
}))

# Median and MAD across samples per protein-gene pair, dataset and metric
plot_data <- results %>%
  group_by(gene, protein, dataset) %>%
  summarise(across(names(metric_labels),
                   list(med = ~median(.x, na.rm = TRUE), mad = ~mad(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            .groups = "drop") %>%
  pivot_longer(-c(gene, protein, dataset),
               names_to = c("metric", ".value"),
               names_pattern = "(.*)_(med|mad)$") %>%
  rename(score = med) %>%
  mutate(metric = factor(recode(metric, !!!metric_labels), levels = metric_labels))

# Protein order: Ward.D2 clustering on medians z-scored per metric (pooling datasets),
# missing dataset-metric combinations imputed with the column median
cluster_matrix <- plot_data %>%
  arrange(protein, dataset, metric) %>%
  group_by(metric) %>%
  mutate(score_scaled = as.numeric(scale(score))) %>%
  ungroup() %>%
  pivot_wider(id_cols = protein, names_from = c(dataset, metric), values_from = score_scaled) %>%
  tibble::column_to_rownames("protein") %>%
  as.matrix()
cluster_matrix <- apply(cluster_matrix, 2, function(x) replace(x, is.na(x), median(x, na.rm = TRUE)))

protein_order <- rownames(cluster_matrix)[hclust(dist(cluster_matrix), method = "ward.D2")$order]
plot_data <- mutate(plot_data, protein = factor(protein, levels = rev(protein_order)))

# Figure 3B: concordance metrics per protein, one colour scale per metric
p <- ggplot(plot_data, aes(x = dataset, y = protein))
metrics <- levels(plot_data$metric)
for (i in seq_along(metrics)) {
  is_last <- i == length(metrics)
  p <- p +
    geom_point(data = filter(plot_data, metric == metrics[i]), aes(colour = score, size = mad), shape = 16) +
    scale_colour_gradient2(name = metrics[i],
                           low = metric_colours[[metrics[i]]][["low"]], mid = "grey95",
                           high = metric_colours[[metrics[i]]][["high"]], midpoint = 0,
                           guide = guide_colourbar(order = i)) +
    scale_size_continuous(name = "MAD", trans = "reverse", guide = if (is_last) "legend" else "none")
  if (!is_last) p <- p + new_scale_colour() + new_scale("size")
}

p <- p +
  facet_grid(cols = vars(metric), scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "grey60", fill = NA, linetype = "dashed"),
        panel.spacing.x = grid::unit(0.15, "cm"),
        panel.grid.major.y = element_line(colour = "grey92"),
        legend.position = "right")

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig3b_concordance_summary.pdf"), p, width = 7, height = 7, units = "in")