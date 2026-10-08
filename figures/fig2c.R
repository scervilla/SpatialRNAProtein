library(dplyr)
library(ggplot2)

path_output <- "figures/figure2/"

datasets <- c("LUAD", "mCRPC")
min_score <- 0.4

# Percentage of expressing cells (intensity > 0) per protein, dataset and cell type
pct_celltype <- bind_rows(lapply(datasets, function(dataset) {
  # Raw protein intensities of all cells (columns named <sample>_<cell>)
  protein <- readRDS(paste0("output/protein_expression/prot_", dataset, ".rds"))
  rownames(protein) <- sub("\\.1$", "", rownames(protein))
  cell_sample <- sub("_.*", "", colnames(protein))
  
  # Cell type labels, low-confidence and unlabelled cells as Unknown
  labels <- bind_rows(lapply(unique(cell_sample), function(sample) {
    predictions <- read.table(paste0("output/cell_annotation/", dataset, "/", sample, "_celltype.txt"))
    data.frame(cell = paste0(sample, "_", rownames(predictions)),
               cell_type = if_else(predictions$prediction.score.max > min_score, predictions$predicted.id, "Unknown"))
  }))
  cell_type <- labels$cell_type[match(colnames(protein), labels$cell)]
  cell_type[is.na(cell_type)] <- "Unknown"
  
  # asinh(x / 10) > 0 is equivalent to x > 0
  positive <- t(as.matrix(protein) > 0) * 1
  pct <- rowsum(positive, cell_type) / as.vector(table(cell_type)) * 100
  
  as.data.frame(as.table(pct)) %>%
    rename(cell_type = Var1, protein = Var2, pct_expr = Freq) %>%
    mutate(dataset = dataset,
           pct_all = (colMeans(positive) * 100)[as.character(protein)])
}))

# Tau index across cell types: 0 = broadly expressed, 1 = restricted to one cell type
tau <- pct_celltype %>%
  group_by(protein, dataset) %>%
  summarise(tau = sum(1 - pct_expr / max(pct_expr)) / (n() - 1),
            pct_expr = first(pct_all),
            .groups = "drop")

# Proteins with tau < 0.1 in any dataset first, then by increasing median tau
protein_order <- tau %>%
  group_by(protein) %>%
  summarise(has_low = any(tau < 0.1, na.rm = TRUE),
            med_tau = median(tau, na.rm = TRUE)) %>%
  arrange(desc(has_low), med_tau) %>%
  pull(protein) %>%
  as.character()

# Figure 2C: tau index and percentage of expressing cells per protein
p <- tau %>%
  mutate(protein = factor(protein, levels = protein_order)) %>%
  ggplot(aes(x = protein, y = dataset, size = pct_expr, colour = tau)) +
  geom_point() +
  scale_color_gradient(low = "gray80", high = "red") +
  labs(x = "Protein", y = "Dataset",
       colour = "Tau Index across\ncell types", size = "% expression") +
  theme_light(base_size = 12) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = "top")

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig2c_tau_dotplot.pdf"), p, width = 7, height = 4, units = "in")