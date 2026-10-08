library(Seurat)
library(ggplot2)

path_output <- "figures/figure2/"

sample <- "TMA20"
marker <- "CD8A"
cofactor <- 10

# Zoomed region (µm, Seurat plot coordinates)
x <- 16938
y <- 2864
half_width <- 50
half_height <- 25
x_range <- c(x - half_width, x + half_width)
y_range <- c(y - half_height, y + half_height)

object <- readRDS(paste0("objects/LUAD/", sample, ".rds"))

# Protein features share names with genes and are suffixed ".1" by Seurat
feature <- paste0(marker, ".1")
stopifnot(feature %in% rownames(object[["Protein"]]))
object$intensity <- asinh(object[["Protein"]]$counts[feature, ] / cofactor)

object[["zoom"]] <- Crop(object[["fov"]], x = x_range, y = y_range, coords = "plot")
DefaultBoundary(object[["zoom"]]) <- "segmentation"

# Figure 2B: asinh-normalized CD8A intensity per cell segmentation
p <- ImageFeaturePlot(object, "intensity", fov = "zoom", border.color = "#FEFFB8", border.size = 1) +
  annotate("rect", xmin = x_range[1], xmax = x_range[2], ymin = y_range[1], ymax = y_range[2],
           colour = "yellow", fill = NA)

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig2b_CD8A_zoom.pdf"), p, width = 5, height = 3, units = "in")