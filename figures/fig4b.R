library(dplyr)
library(ggplot2)

path_output <- "figures/figure4/"

metric_labels <- c(mean_expression = "Mean expression",
                   W0 = "Area",
                   W1 = "Perimeter",
                   W2 = "Euler characteristic",
                   beta = "Anisotropy")

# Kendall's tau-b per dataset, pair and metric (from rank_concordance.R)
kendall <- read.table("output/rank_concordance/kendall.txt") %>%
  filter(metric %in% names(metric_labels)) %>%
  mutate(metric_name = factor(metric_labels[metric], levels = metric_labels),
         metric_type = if_else(metric == "mean_expression", "Pseudobulk", "Spatial"))

# Figure 4B: RNA-protein rank agreement across samples for pseudobulk and spatial metrics
p <- ggplot(kendall, aes(x = metric_name, y = kendall, fill = metric_type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  scale_fill_manual(values = c("Pseudobulk" = "grey40", "Spatial" = "grey90")) +
  labs(x = "", y = "Kendall's tau-b", fill = "Metric type") +
  theme_classic(base_size = 28) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0),
        legend.position = "inside",
        legend.position.inside = c(0.02, 0.02),
        legend.justification = c(0, 0))

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig4b_kendall_boxplot.pdf"), p, width = 6, height = 8, units = "in")
