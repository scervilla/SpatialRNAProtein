library(Seurat)
library(dplyr)
library(ggplot2)

path_output <- "figures/figure3/"

sample <- "TMA20"
marker <- "CD4"
cell_types <- c("T cell CD4", "Macrophages")
min_counts <- 3
min_score <- 0.4
cofactor <- 10

# Zoomed region (µm, Seurat plot coordinates)
x_range <- c(17000, 17400)
y_range <- c(2550, 2700)

object <- readRDS(paste0("objects/LUAD/", sample, ".rds"))
object <- subset(object, subset = nCount_Xenium > min_counts)
object <- NormalizeData(object, assay = "Xenium")

# CD4 protein: background correction and asinh normalization
# Protein features share names with genes and are suffixed ".1" by Seurat
feature <- paste0(marker, ".1")
thresholds <- read.table("utils/thresholds_protein_LUAD.txt")
threshold <- thresholds$threshold[thresholds$protein == marker]
stopifnot(length(threshold) == 1, feature %in% rownames(object[["Protein"]]))
object$protein <- asinh(pmax(object[["Protein"]]$counts[feature, ] - threshold, 0) / cofactor)

# CD4 T cells and macrophages; other and low-confidence cells as "na"
predictions <- read.table(paste0("output/cell_annotation/LUAD/", sample, "_celltype.txt"))
stopifnot(identical(rownames(predictions), colnames(object)))
object$cell_type <- if_else(predictions$predicted.id %in% cell_types & predictions$prediction.score.max >= min_score,
                            predictions$predicted.id, "na")

object[["zoom"]] <- Crop(object[["fov"]], x = x_range, y = y_range, coords = "plot")
DefaultBoundary(object[["zoom"]]) <- "segmentation"
zoom_box <- annotate("rect", xmin = x_range[1], xmax = x_range[2], ymin = y_range[1], ymax = y_range[2],
                     colour = "yellow", fill = NA)

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)

# Figure 3D: CD4 T cells and macrophages
p <- ImageDimPlot(object, group.by = "cell_type", fov = "zoom", border.color = "#FEFFB8", border.size = 0.1, axes = TRUE) +
  scale_fill_manual(values = c("na" = "gray20", "Macrophages" = "yellow2", "T cell CD4" = "purple3")) +
  zoom_box +
  NoLegend()
ggsave(paste0(path_output, "fig3d_celltype.pdf"), p, width = 20, height = 8, units = "in")

# Figure 3D: CD4 RNA (log-normalized) and protein (background-corrected asinh); protein also with legend
feature_plots <- data.frame(feature = c(marker, "protein", "protein"),
                            name = c("CD4_rna", "CD4_protein", "CD4_protein_legend"),
                            legend = c(FALSE, FALSE, TRUE))
for (i in seq_len(nrow(feature_plots))) {
  p <- ImageFeaturePlot(object, features = feature_plots$feature[i], fov = "zoom",
                        border.color = "#FEFFB8", border.size = 0.1, axes = TRUE) +
    scale_fill_gradient(low = "gray20", high = "purple") +
    zoom_box
  if (!feature_plots$legend[i]) p <- p + NoLegend()
  ggsave(paste0(path_output, "fig3d_", feature_plots$name[i], ".pdf"), p, width = 20, height = 8, units = "in")
}