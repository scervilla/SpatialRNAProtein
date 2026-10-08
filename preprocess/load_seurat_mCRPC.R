library(progressr)
library(Seurat)
library(dplyr)

ReadXenium <- function (data.dir, outs = c("matrix"), type = "centroids", 
                        mols.qv.threshold = 20) {
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"), 
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix"), 
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      pmtx <- progressor()
      pmtx(message = "Reading counts matrix", class = "sticky", 
           amount = 0)
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir, 
                                                              "cell_feature_matrix/")))
      pmtx(type = "finish")
      matrix
    }, centroids = {
      pcents <- progressor()
      pcents(message = "Loading cell centroids", class = "sticky", 
             amount = 0)
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid, 
                                     y = cell_info$y_centroid, cell = cell_info$cell_id, 
                                     stringsAsFactors = FALSE)
      pcents(type = "finish")
      cell_centroid_df
    }, segmentations = {
      psegs <- progressor()
      psegs(message = "Loading cell segmentations", class = "sticky", 
            amount = 0)
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir, 
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      psegs(type = "finish")
      cell_boundaries_df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}

LoadXenium <- function (data.dir, fov = "fov", assay = "Xenium") 
{
  
  data <- ReadXenium(data.dir = data.dir, type = c("centroids", 
                                                   "segmentations"), )
  segmentations.data <- list(centroids = CreateCentroids(data$centroids), 
                             segmentation = CreateSegmentation(data$segmentations))
  coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                            "centroids"), assay = assay)
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], 
                                   assay = assay)
  if ("Protein Expression" %in% names(data$matrix)) {
    xenium.obj[["Protein"]] <- CreateAssayObject(counts = data$matrix[["Protein Expression"]]) 
  }
  if ("Blank Codeword" %in% names(data$matrix)) 
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
  else xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}


file_names <- c(
  "output-XETG00289__0102369__FISPLAT18__20260604__144148",
  "output-XETG00289__0102369__FISPLAT29__20260604__144148",
  "output-XETG00289__0102369__FISPLAT5__20260604__144148",
  "output-XETG00289__0102369__FISPLAT7__20260604__144148",
  "output-XETG00289__0102369__FISPLAT8__20260604__144148",
  "output-XETG00289__0102369__FISPLAT9__20260604__144148",
  "output-XETG00289__0102372__FISPLAT2__20260604__144148",
  "output-XETG00289__0102372__FISPLAT22__20260604__144148",
  "output-XETG00289__0102372__FISPLAT26__20260604__144148",
  "output-XETG00289__0102372__FISPLAT28__20260604__144148",
  "output-XETG00289__0102372__FISPLAT30__20260604__144148"
)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

path_raw <- "raw/mCRPC/"
path_objects <- "objects/mCRPC/"

names(raw_data) <- stringr::str_split_i(raw_data, "__", 3)
stopifnot(sample %in% names(raw_data))

object <- LoadXenium(paste0(path_raw, raw_data[sample], "/"))

dir.create(path_objects, recursive = TRUE, showWarnings = FALSE)
saveRDS(object, paste0(path_objects, sample, ".rds"))

