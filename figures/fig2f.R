library(Seurat)
library(dplyr)
library(ggplot2)

path_output <- "figures/figure2/"

sample <- "TMA20"
marker <- "CD68"
cofactor <- 10

object <- readRDS(paste0("objects/LUAD/", sample, ".rds"))
thresholds <- read.table("utils/thresholds_protein_LUAD.txt")
threshold <- thresholds$threshold[thresholds$protein == marker]
stopifnot(length(threshold) == 1)

# Protein features share names with genes and are suffixed ".1" by Seurat
feature <- paste0(marker, ".1")
stopifnot(feature %in% rownames(object[["Protein"]]))
intensity <- object[["Protein"]]$counts[feature, ]

# asinh-normalized intensity before and after background correction
object$original <- asinh(intensity / cofactor)
object$corrected <- asinh(pmax(intensity - threshold, 0) / cofactor)

# Figure 2F (spatial): CD68 intensity per cell segmentation before and after correction
DefaultBoundary(object[["fov"]]) <- "segmentation"
p_spatial <- ImageFeaturePlot(object, c("original", "corrected"), fov = "fov", border.size = 0) &
  scale_fill_gradient(low = "gray20", high = "green") &
  NoLegend()

# Figure 2F (histogram): intensity distribution before and after correction
plot_data <- data.frame(variable = factor(rep(c("Original", "Corrected"), each = ncol(object)),
                                          levels = c("Original", "Corrected")),
                        value = c(object$original, object$corrected))
zero_stats <- plot_data %>%
  group_by(variable) %>%
  summarise(pct_zero = mean(value == 0) * 100)

p_histogram <- ggplot(plot_data, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "grey40", colour = "white") +
  geom_text(data = zero_stats, aes(x = 0.25, y = Inf, label = sprintf("Zeros: %.1f%%", pct_zero)),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.3, size = 3.5) +
  facet_wrap(~variable, scales = "free", ncol = 2) +
  scale_y_sqrt() +
  labs(x = "Normalized intensity", y = "Density (sqrt scale)", title = paste(marker, "in", sample)) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig2f_CD68_spatial.pdf"), p_spatial, width = 6, height = 3, units = "in")
ggsave(paste0(path_output, "fig2f_CD68_histogram.pdf"), p_histogram, width = 5.5, height = 2.75, units = "in")