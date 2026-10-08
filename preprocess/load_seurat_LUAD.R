library(Seurat)
library(dplyr)
library(arrow)

args = commandArgs(trailingOnly=TRUE)
slide_id <- args[1]

path_raw <- "raw/LUAD/"
path_objects <- "objects/LUAD/"

# Load Xenium gene and protein matrices with cell centroids and segmentations
load_xenium_protein <- function(data_dir, fov = "fov", assay = "Xenium") {
  matrix <- Read10X(data.dir = file.path(data_dir, "cell_feature_matrix/"))
  stopifnot(all(c("Gene Expression", "Protein Expression") %in% names(matrix)))

  cell_info <- data.table::fread(file.path(data_dir, "cells.csv.gz"))
  cell_boundaries <- data.table::fread(file.path(data_dir, "cell_boundaries.csv.gz"))

  centroids <- data.frame(x = cell_info$x_centroid, y = cell_info$y_centroid, cell = cell_info$cell_id)
  segmentations <- data.frame(cell = cell_boundaries$cell_id, x = cell_boundaries$vertex_x, y = cell_boundaries$vertex_y)

  coords <- CreateFOV(coords = list(centroids = CreateCentroids(centroids),
                                    segmentation = CreateSegmentation(segmentations)),
                      type = c("segmentation", "centroids"), assay = assay)

  object <- CreateSeuratObject(counts = matrix[["Gene Expression"]], assay = assay)
  object[["Protein"]] <- CreateAssayObject(counts = matrix[["Protein Expression"]])
  blank <- if ("Blank Codeword" %in% names(matrix)) "Blank Codeword" else "Unassigned Codeword"
  object[["BlankCodeword"]] <- CreateAssayObject(counts = matrix[[blank]])
  object[["ControlCodeword"]] <- CreateAssayObject(counts = matrix[["Negative Control Codeword"]])
  object[["ControlProbe"]] <- CreateAssayObject(counts = matrix[["Negative Control Probe"]])
  object[[fov]] <- coords
  object
}

object <- load_xenium_protein(paste0(path_raw, slide_id, "/"))

# Add TMA core assignment (from cells_TMA.R)
cells_tma <- read_parquet(paste0(path_objects, "cells_TMA.parquet")) %>%
  filter(slide == slide_id) %>%
  select(cell_id, cell_area:TMA)
stopifnot(nrow(cells_tma) > 0, all(cells_tma$cell_id %in% colnames(object)))

object <- AddMetaData(object, tibble::column_to_rownames(as.data.frame(cells_tma), "cell_id"))
saveRDS(object, paste0(path_objects, slide_id, ".rds"))

# Split by TMA core, cells outside any core are dropped
object_list <- SplitObject(subset(object, cells = cells_tma$cell_id), split.by = "TMA")

dir.create(paste0(path_objects), recursive = TRUE, showWarnings = FALSE)
for (tma in names(object_list)) {
  saveRDS(object_list[[tma]], paste0(path_objects, tma, ".rds"))
}
