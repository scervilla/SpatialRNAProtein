library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

path_objects <- "objects/LUAD/"
path_panel <- "utils/gene_panels/"
path_reference <- "reference/"
path_output <- "output/cell_annotation/LUAD/"

min_counts_query <- 3
min_counts_reference <- 10
n_dims <- 30

excluded_types <- c("transitional club/AT2", "other", "Alveolar cell type 1", "Club")
coarse_labels <- c("Alveolar cell type 2" = "Epithelial/Tumor",
                   "Tumor cells"          = "Epithelial/Tumor",
                   "Ciliated"             = "Epithelial/Tumor",
                   "Endothelial cell"     = "Endothelial cells",
                   "Stromal"              = "Fibroblasts",
                   "Plasma cell"          = "Plasma cells",
                   "B cell"               = "B cells",
                   "T cell CD4"           = "T cell CD4",
                   "T cell CD8"           = "T cell CD8",
                   "T cell regulatory"    = "T cell regulatory",
                   "NK cell"              = "NK cells",
                   "Macrophage"           = "Macrophages",
                   "Macrophage alveolar"  = "Macrophages",
                   "Monocyte"             = "Monocytes",
                   "Neutrophils"          = "Neutrophils",
                   "Mast cell"            = "Mast cells",
                   "cDC1"                 = "DC",
                   "cDC2"                 = "DC",
                   "DC mature"            = "DC",
                   "pDC"                  = "pDC")

# Reference: subset to Xenium panel genes and rename Ensembl IDs to gene symbols
panel <- read.csv(paste0(path_panel, "Xenium_hIO_v1_metadata.csv"))
reference <- readRDS(paste0(path_reference, "lung_sc_subset.rds"))
panel <- panel[panel$Ensemble.ID %in% rownames(reference), ]
reference <- reference[panel$Ensemble.ID, ]
rownames(reference) <- panel$Gene

# Remove low-quality cells and cell types not resolved by the panel
reference$rna_count <- colSums(reference[["RNA"]]$counts)
reference <- subset(reference, subset = rna_count >= min_counts_reference & !cell_type_major %in% excluded_types)

# Collapse to coarse labels, pericytes from the fine annotation
reference$celltype_coarse <- unname(coarse_labels[as.character(reference$cell_type_major)])
reference$celltype_coarse[reference$cell_type == "pericyte"] <- "Pericytes"
stopifnot(!anyNA(reference$celltype_coarse))

reference <- NormalizeData(reference) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = n_dims, verbose = FALSE)

# Query: TMA core
object <- readRDS(paste0(path_objects, sample, ".rds"))
object <- subset(object, subset = nCount_Xenium > min_counts_query)
DefaultAssay(object) <- "Xenium"
object <- NormalizeData(object) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = n_dims, verbose = FALSE)

# Label transfer
options(future.globals.maxSize = 1.5 * 1024^3)
common_genes <- intersect(rownames(reference), rownames(object))
anchors <- FindTransferAnchors(reference = reference,
                               query = object,
                               query.assay = "Xenium",
                               features = common_genes,
                               reference.reduction = "pca",
                               dims = 1:n_dims)
predictions <- TransferData(anchorset = anchors,
                            refdata = reference$celltype_coarse,
                            dims = 1:n_dims)

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
write.table(predictions, paste0(path_output, sample, "_celltype.txt"))
