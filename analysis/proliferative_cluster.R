library(Seurat)
library(dplyr)
library(readxl)
library(Matrix)
library(scPearsonPCA)

path_objects <- "objects/LUAD/"
path_output <- "output/proliferative/LUAD/"

patient <- 1
min_counts <- 3
n_dims <- 30
clip_sd <- 10
cofactor <- 10
proliferative_cluster <- "10"

set.seed(1)

# TMA cores of the selected patient
samples <- readLines("utils/samples_LUAD.txt")
cores <- read_xlsx("utils/TMA_metadata.xlsx") %>%
  filter(tma_id %in% samples, patient_id == patient) %>%
  pull(tma_id)
stopifnot(length(cores) > 0)

objects <- lapply(cores, function(core) readRDS(paste0(path_objects, core, ".rds")))
object <- merge(objects[[1]], objects[-1], add.cell.ids = cores)
object <- JoinLayers(object)
object <- subset(object, subset = nCount_Xenium > min_counts)
DefaultAssay(object) <- "Xenium"

# Highly variable genes with non-zero counts
object <- FindVariableFeatures(object)
counts <- object[["Xenium"]]$counts
hvgs <- VariableFeatures(object)
hvgs <- hvgs[rowSums(counts[hvgs, ]) > 0]

# Sparse quasi-Poisson PCA on Pearson residuals, clipped at 10 SD, centred and scaled
pca <- sparse_quasipoisson_pca_seurat(counts[hvgs, ],
                                      totalcounts = colSums(counts),
                                      grate = gene_frequency(counts)[hvgs],
                                      scale.max = clip_sd,
                                      do.scale = TRUE,
                                      do.center = TRUE)
object[["pca"]] <- pca$reduction.data

# UMAP, shared nearest-neighbour graph and Louvain clustering
object <- NormalizeData(object) %>%
  RunUMAP(dims = 1:n_dims) %>%
  FindNeighbors(dims = 1:n_dims) %>%
  FindClusters()

# Cluster markers and differential expression of the proliferative cluster against all other cells
markers_clusters <- FindAllMarkers(object, only.pos = TRUE)
mki67_cluster <- tapply(FetchData(object, "MKI67")[, 1], Idents(object), mean)
stopifnot(names(which.max(mki67_cluster)) == proliferative_cluster)
markers_proliferative <- FindMarkers(object, ident.1 = proliferative_cluster)
markers_proliferative$gene <- rownames(markers_proliferative)

# Protein background correction and asinh normalization (for visualization)
protein <- as.matrix(object[["Protein"]]$counts)
thresholds <- read.table("utils/thresholds_protein_LUAD.txt")
idx <- match(thresholds$protein, sub("\\.1$", "", rownames(protein)))
stopifnot(!anyNA(idx))
protein[idx, ] <- pmax(protein[idx, ] - thresholds$threshold, 0)
object[["Protein"]]$data <- as(asinh(protein / cofactor), "dgCMatrix")

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
write.table(markers_clusters, paste0(path_output, "markers_clusters_patient", patient, ".txt"))
write.table(markers_proliferative, paste0(path_output, "DEG_proliferative_patient", patient, ".txt"))
saveRDS(object, paste0(path_objects, "patient", patient, ".rds"))