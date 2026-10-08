library(Seurat)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

path_objects <- "objects/mCRPC/"
path_panel <- "utils/gene_panels/"
path_reference <- "reference/"
path_output <- "output/cell_annotation/mCRPC/"

min_counts_query <- 3
min_counts_reference <- 10
n_dims <- 30

excluded_types <- c("Cycling T")
coarse_labels <- c("Epitheial_Luminal"   = "Epithelial/Tumor",
                   "Epitheial_Basal"     = "Epithelial/Tumor",
                   "Epitheial_Club"      = "Epithelial/Tumor",
                   "Epitheial_Hillock"   = "Epithelial/Tumor",
                   "Tumor"               = "Epithelial/Tumor",
                   "Endothelial cells-1" = "Endothelial cells",
                   "Endothelial cells-2" = "Endothelial cells",
                   "Fibroblasts"         = "Fibroblasts",
                   "Pericytes-1"         = "Pericytes",
                   "Pericytes-2"         = "Pericytes",
                   "Plasma cells"        = "Plasma cells",
                   "B cells"             = "B cells",
                   "Naive Th"            = "T cell CD4",
                   "Th1"                 = "T cell CD4",
                   "Th17"                = "T cell CD4",
                   "CTL-1"               = "T cell CD8",
                   "CTL-2"               = "T cell CD8",
                   "Treg"                = "T cell regulatory",
                   "TNK"                 = "NKT cells",
                   "NK"                  = "NK cells",
                   "Macrophage1"         = "Macrophages",
                   "Macrophage2"         = "Macrophages",
                   "Macrophage3"         = "Macrophages",
                   "Mono1"               = "Monocytes",
                   "Mono2"               = "Monocytes",
                   "Mono3"               = "Monocytes",
                   "Mast cells"          = "Mast cells",
                   "mDC"                 = "DC",
                   "PDC"                 = "pDC")

# Reference: RNA assay subset to Xenium panel genes
panel <- read.csv(paste0(path_panel, "Xenium_hIO_v1_metadata.csv"))
reference <- readRDS(paste0(path_reference, "GSE181294_scRNA_reference_light.rds"))
DefaultAssay(reference) <- "RNA"
reference[["SCT"]] <- NULL
reference <- JoinLayers(reference[intersect(panel$Gene, rownames(reference)), ])

# Remove low-quality cells and cell types not resolved by the panel
reference$rna_count <- colSums(reference[["RNA"]]$counts)
reference <- subset(reference, subset = rna_count >= min_counts_reference & !celltype %in% excluded_types)

# Collapse to coarse labels
reference$celltype_coarse <- unname(coarse_labels[as.character(reference$celltype)])
stopifnot(!anyNA(reference$celltype_coarse))

reference <- NormalizeData(reference) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = n_dims, verbose = FALSE)

# Query: mCRPC sample
object <- readRDS(paste0(path_objects, sample, ".rds"))
object <- subset(object, subset = nCount_Xenium > min_counts_query)
DefaultAssay(object) <- "Xenium"
object <- NormalizeData(object) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = n_dims, verbose = FALSE)

# Label transfer
options(future.globals.maxSize = 5.2 * 1024^3)
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
