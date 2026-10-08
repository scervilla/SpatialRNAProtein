library(ggplot2)

path_output <- "figures/figure2/"

cofactor <- 10
marker <- "CD8A"

# Raw protein intensities of all LUAD cells (columns named <sample>_<cell>)
protein <- readRDS("output/protein_expression/prot_LUAD.rds")
samples <- readLines("utils/samples_LUAD.txt")

# Protein features share names with genes and are suffixed ".1" by Seurat
rownames(protein) <- sub("\\.1$", "", rownames(protein))
stopifnot(marker %in% rownames(protein))

protein <- protein[, sub("_.*", "", colnames(protein)) %in% samples]
stopifnot(ncol(protein) > 0)

# Figure 2A: distribution of asinh-normalized CD8A intensity
plot_data <- data.frame(intensity = asinh(protein[marker, ] / cofactor))

p <- ggplot(plot_data, aes(x = intensity)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "gray40", colour = "white") +
  scale_y_sqrt(expand = expansion(mult = c(0.01, 0.05))) +
  labs(x = "Normalized intensity", y = "Density (sqrt scale)", title = marker) +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig2a_CD8A_histogram.pdf"), p, width = 2.75, height = 3.5, units = "in")
