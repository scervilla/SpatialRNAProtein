library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

path_output <- "figures/figure2/"

datasets <- c("LUAD", "mCRPC")
min_score <- 0.4
min_diff_label <- 0.1

# Percentage of expressing cells (intensity > 0) per protein, dataset and cell type,
# before and after background correction
pct_celltype <- bind_rows(lapply(datasets, function(dataset) {
  # Raw protein intensities of all cells (columns named <sample>_<cell>)
  protein <- as.matrix(readRDS(paste0("output/protein_expression/prot_", dataset, ".rds")))
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
  
  # Background correction (proteins without a threshold are left unchanged)
  thresholds <- read.table(paste0("utils/thresholds_protein_", dataset, ".txt"))
  stopifnot(all(thresholds$protein %in% rownames(protein)))
  corrected <- protein
  corrected[thresholds$protein, ] <- pmax(corrected[thresholds$protein, ] - thresholds$threshold, 0)
  
  # asinh(x / 10) > 0 is equivalent to x > 0
  bind_rows(lapply(list(original = protein, corrected = corrected), function(intensity) {
    pct <- rowsum(t(intensity > 0) * 1, cell_type) / as.vector(table(cell_type)) * 100
    as.data.frame(as.table(pct)) %>%
      rename(cell_type = Var1, protein = Var2, pct_expr = Freq)
  }), .id = "expression") %>%
    mutate(dataset = dataset)
}))

# Tau index across cell types: 0 = broadly expressed, 1 = restricted to one cell type
tau <- pct_celltype %>%
  group_by(protein, dataset, expression) %>%
  summarise(tau = sum(1 - pct_expr / max(pct_expr)) / (n() - 1), .groups = "drop") %>%
  pivot_wider(names_from = expression, values_from = tau) %>%
  mutate(diff = abs(corrected - original))

# Figure 2G: tau index before and after background correction
p <- ggplot(tau, aes(x = original, y = corrected)) +
  geom_point() +
  geom_label_repel(data = filter(tau, diff > min_diff_label), aes(label = protein), size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~dataset) +
  labs(x = "Tau Index in Original expression",
       y = "Tau Index in Corrected expression") +
  theme_classic(base_size = 12) +
  theme(aspect.ratio = 1)

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig2g_tau_correction.pdf"), p, width = 7, height = 4, units = "in")