library(dplyr)
library(tidyr)
library(ggplot2)

path_output <- "figures/figure3/"

datasets <- c("LUAD", "mCRPC")
min_cells <- 10
min_rank_diff <- 4

# Percentage of RNA- and protein-positive cells per cell type and sample (from expression_celltype.R)
results <- bind_rows(lapply(datasets, function(dataset) {
  samples <- readLines(paste0("utils/samples_", dataset, ".txt"))
  bind_rows(lapply(samples, function(sample) {
    read.table(paste0("output/expression_celltype/", dataset, "/", sample, "_summary.txt")) %>%
      mutate(sample = sample)
  })) %>%
    mutate(dataset = dataset)
}))

# Rank difference of each cell type between RNA and protein positivity, per sample and marker;
# summarised per cell type and marker as the fraction of samples with |rank diff| > 4 and the median rank diff
plot_data <- results %>%
  filter(n > min_cells, cell_type != "Unknown") %>%
  group_by(sample, protein, dataset) %>%
  mutate(rank_diff = rank(pct_rna_pos) - rank(pct_prot_pos)) %>%
  group_by(cell_type, protein) %>%
  summarise(frac_samples = sum(abs(rank_diff) > min_rank_diff) / n_distinct(sample),
            median_rank_diff = median(rank_diff),
            .groups = "drop")

# Cell type and marker order: hierarchical clustering on the median rank difference (missing pairs = 0)
mat <- plot_data %>%
  select(cell_type, protein, median_rank_diff) %>%
  pivot_wider(names_from = protein, values_from = median_rank_diff, values_fill = 0) %>%
  tibble::column_to_rownames("cell_type") %>%
  as.matrix()
cell_order <- rownames(mat)[hclust(dist(mat))$order]
protein_order <- colnames(mat)[hclust(dist(t(mat)))$order]

# Figure 3C: cell type ranking discordance between RNA and protein per marker
p <- plot_data %>%
  mutate(cell_type = factor(cell_type, levels = cell_order),
         protein = factor(protein, levels = protein_order)) %>%
  ggplot(aes(x = protein, y = cell_type, size = frac_samples, colour = median_rank_diff)) +
  geom_point() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "#B2182B") +
  labs(x = "", y = "Cell Type", size = "% samples", colour = "Rank diff") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        panel.border = element_rect(colour = "grey60", fill = NA, linetype = "dashed"),
        panel.spacing.x = grid::unit(0.15, "cm"),
        panel.grid.major.y = element_line(colour = "grey92"),
        panel.grid.major.x = element_line(colour = "grey92"),
        legend.position = "top")

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig3c_celltype_rank_diff.pdf"), p, width = 8, height = 5.5, units = "in")