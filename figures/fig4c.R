library(dplyr)
library(readxl)
library(ggplot2)

path_output <- "figures/figure4/"

pair <- "PD-L1"
dataset <- "LUAD"
min_samples <- 3
metric_labels <- c(mean_expression = "Mean expression",
                   W0 = "Area",
                   W1 = "Perimeter",
                   W2 = "Euler characteristic",
                   beta = "Anisotropy")

# Paired RNA and protein Minkowski values per sample (from rank_concordance.R), with core region
regions <- read_xlsx("utils/TMA_metadata.xlsx") %>% select(sample = tma_id, region)
paired <- read.table("output/rank_concordance/minkowski_paired.txt") %>%
  filter(pair == !!pair, dataset == !!dataset) %>%
  left_join(regions, by = "sample")

# Figure 4C: RNA vs protein per sample for each metric, all cores and tumour-region cores
subsets <- list(all = paired, tumour = filter(paired, region == "T"))
dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
for (subset_name in names(subsets)) {
  for (metric in names(metric_labels)) {
    plot_data <- filter(subsets[[subset_name]], metric == !!metric)
    n <- sum(complete.cases(plot_data$rna, plot_data$protein))
    tau <- if (n >= min_samples) round(cor(plot_data$rna, plot_data$protein, method = "kendall", use = "complete.obs"), 3) else NA

    p <- ggplot(plot_data, aes(x = rna, y = protein)) +
      geom_point() +
      scale_x_log10() +
      scale_y_log10() +
      labs(x = "RNA (log10)", y = "Protein (log10)",
           title = paste0(metric_labels[metric], " ", pair, "\n(τ = ", tau, ")")) +
      theme_classic(base_size = 22) +
      theme(aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5, size = 25))

    suffix <- if (subset_name == "tumour") "_tumour" else ""
    ggsave(paste0(path_output, "fig4c_PDL1_", metric, suffix, ".pdf"), p, width = 6, height = 6, units = "in")
  }
}
