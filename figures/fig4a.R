library(Seurat)
library(dplyr)
library(ggplot2)
library(sf)
library(igraph)

path_output <- "figures/figure4/"

sample <- "TMA20"
min_counts <- 3
min_score <- 0.4
cofactor <- 10
contact_tol <- 3     # µm buffer for cell contact between tumor segmentations
min_comp_size <- 10  # minimum cells per tumor component
union_tol <- 5       # µm buffer to merge touching segmentations into one outline
pairs <- data.frame(gene = c("LAG3", "CD274"), protein = c("LAG-3", "PD-L1"), name = c("LAG3", "PDL1"))

object <- readRDS(paste0("objects/LUAD/", sample, ".rds"))
object <- subset(object, subset = nCount_Xenium > min_counts)
object <- NormalizeData(object, assay = "Xenium")

# Protein background correction and asinh normalization
protein <- as.matrix(object[["Protein"]]$counts)
thresholds <- read.table("utils/thresholds_protein_LUAD.txt")
idx <- match(thresholds$protein, sub("\\.1$", "", rownames(protein)))
stopifnot(!anyNA(idx))
protein[idx, ] <- pmax(protein[idx, ] - thresholds$threshold, 0)
object[["Protein"]]$data <- as(asinh(protein / cofactor), "dgCMatrix")

# tumor cells (low-confidence labels excluded)
predictions <- read.table(paste0("output/cell_annotation/LUAD/", sample, "_celltype.txt"))
stopifnot(identical(rownames(predictions), colnames(object)))
tumor_cells <- colnames(object)[predictions$predicted.id == "Epithelial/Tumor" &
                                   predictions$prediction.score.max >= min_score]

# tumor cell segmentations as closed polygons
segmentation <- GetTissueCoordinates(object[["fov"]], which = "segmentation") %>%
  filter(cell %in% tumor_cells)
cells_sf <- st_sf(cell_id = unique(segmentation$cell),
                  geometry = st_sfc(lapply(split(segmentation[, c("x", "y")], segmentation$cell)[unique(segmentation$cell)],
                                           function(xy) st_polygon(list(as.matrix(rbind(xy, xy[1,])))))))

# Connected tumor components: segmentations within 3 µm of each other
touching <- st_intersects(st_buffer(cells_sf, contact_tol), cells_sf)
edges <- data.frame(from = rep(seq_along(touching), lengths(touching)), to = unlist(touching)) %>%
  filter(from < to)
graph <- graph_from_data_frame(data.frame(from = cells_sf$cell_id[edges$from], to = cells_sf$cell_id[edges$to]),
                               directed = FALSE, vertices = data.frame(name = cells_sf$cell_id))
membership <- components(graph)$membership
large_components <- names(which(table(membership) >= min_comp_size))

# Outline per component: union of segmentations, holes removed
outlines <- do.call(c, lapply(large_components, function(comp) {
  cells_sf %>%
    filter(cell_id %in% names(membership)[membership == as.integer(comp)]) %>%
    st_buffer(union_tol) %>%
    st_union() %>%
    st_make_valid() %>%
    st_buffer(-union_tol) %>%
    st_make_valid() %>%
    st_cast("POLYGON", warn = FALSE) %>%
    lapply(function(g) st_polygon(list(g[[1]]))) %>%
    st_sfc()
}))

# Outline vertices in Seurat plot coordinates (x and y swapped)
outline_df <- as.data.frame(st_coordinates(outlines)) %>%
  transmute(x = Y, y = X, polygon = L2)

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)

# Figure 4A: RNA and protein expression with tumor component outlines
for (i in seq_len(nrow(pairs))) {
  p <- ImageFeaturePlot(object, features = c(pairs$gene[i], pairs$protein[i]), fov = "fov",
                        border.color = NA, border.size = 0) &
    scale_fill_viridis_c(option = "G", trans = "sqrt") &
    geom_polygon(data = outline_df, aes(x = x, y = y, group = polygon), inherit.aes = FALSE,
                 fill = NA, colour = "yellow", linewidth = 0.45) &
    NoLegend()
  ggsave(paste0(path_output, "fig4a_", pairs$name[i], ".png"), p, width = 6, height = 3, units = "in", dpi = 300)
}