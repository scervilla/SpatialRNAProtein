library(Seurat)
library(dplyr)
library(ggplot2)
library(arrow)

path_output <- "figures/figure3/"

sample <- "TMA20"
slide <- "TMA1-SK11-1"
cofactor <- 10
min_qv <- 20
genes <- c("CD8A", "LAG3")

# Zoomed region (µm, Seurat plot coordinates: plot x = y_location, plot y = x_location)
x_range <- c(16850.19, 16912.02)
y_range <- c(2690.037, 2734.237)

object <- readRDS(paste0("objects/LUAD/", sample, ".rds"))
object <- NormalizeData(object, assay = "Xenium")

# Protein background correction and asinh normalization
protein <- as.matrix(object[["Protein"]]$counts)
thresholds <- read.table("utils/thresholds_protein_LUAD.txt")
idx <- match(thresholds$protein, sub("\\.1$", "", rownames(protein)))
stopifnot(!anyNA(idx))
protein[idx, ] <- pmax(protein[idx, ] - thresholds$threshold, 0)
object[["Protein"]]$data <- as(asinh(protein / cofactor), "dgCMatrix")

object[["zoom"]] <- Crop(object[["fov"]], x = x_range, y = y_range, coords = "plot")
DefaultBoundary(object[["zoom"]]) <- "segmentation"
zoom_box <- annotate("rect", xmin = x_range[1], xmax = x_range[2], ymin = y_range[1], ymax = y_range[2],
                     colour = "yellow", fill = NA)

# High-quality CD8A and LAG3 transcripts inside the zoomed region
transcripts <- read_parquet(paste0("raw/LUAD/", slide, "/transcripts.parquet"), as_data_frame = FALSE) %>%
  filter(y_location >= x_range[1], y_location <= x_range[2],
         x_location >= y_range[1], x_location <= y_range[2],
         qv >= min_qv, feature_name %in% genes) %>%
  collect()

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)

# Figure 3A: CD8A and LAG3 transcripts over cell segmentations
p <- ImageDimPlot(object, group.by = "TMA", fov = "zoom", border.size = 0.5, border.color = "#FEFFB8", axes = FALSE) +
  scale_fill_manual(values = "gray20") +
  geom_point(data = transcripts, aes(x = y_location, y = x_location, colour = feature_name),
             inherit.aes = FALSE, size = 4) +
  scale_colour_manual(values = c("CD8A" = "red", "LAG3" = "yellow")) +
  zoom_box
ggsave(paste0(path_output, "fig3a_transcripts.pdf"), p, width = 9, height = 9, units = "in")

# Figure 3A: RNA and protein expression per cell (RNA log-normalized, protein background-corrected asinh)
feature_plots <- data.frame(feature = c("CD8A", "CD8A.1", "LAG3", "LAG-3"),
                            colour = c("red", "red", "yellow", "yellow"),
                            name = c("CD8A_rna", "CD8A_protein", "LAG3_rna", "LAG3_protein"))
for (i in seq_len(nrow(feature_plots))) {
  p <- ImageFeaturePlot(object, features = feature_plots$feature[i], fov = "zoom",
                        border.size = 0.5, border.color = "#FEFFB8", axes = FALSE) +
    scale_fill_gradient(low = "gray20", high = feature_plots$colour[i]) +
    zoom_box
  ggsave(paste0(path_output, "fig3a_", feature_plots$name[i], ".pdf"), p, width = 9, height = 9, units = "in")
}

# Figure 3A: CD8A and LAG3 positive cells (> 0) per modality
positivity <- list(rna = object[["Xenium"]]$data[c("CD8A", "LAG3"), ],
                   protein = object[["Protein"]]$data[c("CD8A.1", "LAG-3"), ])
for (modality in names(positivity)) {
  cd8a_pos <- positivity[[modality]][1, ] > 0
  lag3_pos <- positivity[[modality]][2, ] > 0
  object$positive <- case_when(cd8a_pos & lag3_pos ~ "pos",
                               lag3_pos ~ "lag3",
                               cd8a_pos ~ "cd8a",
                               TRUE ~ "neg")
  p <- ImageDimPlot(object, group.by = "positive", fov = "zoom", border.size = 0.5, border.color = "#FEFFB8", axes = FALSE) +
    scale_fill_manual(values = c("neg" = "gray20", "cd8a" = "red", "lag3" = "yellow", "pos" = "orange")) +
    zoom_box
  ggsave(paste0(path_output, "fig3a_", modality, "_positive.pdf"), p, width = 9, height = 9, units = "in")
}