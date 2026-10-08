library(Seurat)
library(dplyr)
library(readxl)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
dataset <- args[2]

path_objects <- paste0("objects/", dataset, "/")
path_annotation <- paste0("output/cell_annotation/", dataset, "/")
path_output <- paste0("output/expression_celltype/", dataset, "/")

min_counts <- 3
min_score <- 0.4
cofactor <- 10

# Protein-gene pairs (proteins without a panel gene are dropped)
pairs <- read_xlsx("utils/protein_gene.xlsx") %>% filter(!is.na(gene))
thresholds <- read.table(paste0("utils/thresholds_protein_", dataset, ".txt"))

object <- readRDS(paste0(path_objects, sample, ".rds"))
object <- subset(object, subset = nCount_Xenium > min_counts)
object <- NormalizeData(object)

predictions <- read.table(paste0(path_annotation, sample, "_celltype.txt"))
stopifnot(identical(rownames(predictions), colnames(object)))
cell_type <- if_else(predictions$prediction.score.max < min_score, "Unknown", predictions$predicted.id)

# Protein background correction and asinh normalization
# Protein features share names with genes and are suffixed ".1" by Seurat
protein <- as.matrix(object[["Protein"]]$counts)
rownames(protein) <- sub("\\.1$", "", rownames(protein))
stopifnot(all(thresholds$protein %in% rownames(protein)),
          all(pairs$protein %in% rownames(protein)),
          all(pairs$gene %in% rownames(object[["Xenium"]])))

protein[thresholds$protein, ] <- pmax(protein[thresholds$protein, ] - thresholds$threshold, 0)
protein <- asinh(protein / cofactor)
rna <- object[["Xenium"]]$data

# Mean expression and percentage of expressing cells per cell type
summary <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  data.frame(cell_type = cell_type,
             rna = rna[pairs$gene[i], ],
             prot = protein[pairs$protein[i], ]) %>%
    group_by(cell_type) %>%
    summarise(n = n(),
              mean_rna = mean(rna),
              pct_rna_pos = mean(rna > 0) * 100,
              mean_prot = mean(prot),
              pct_prot_pos = mean(prot > 0) * 100) %>%
    mutate(gene = pairs$gene[i], protein = pairs$protein[i])
}))

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
write.table(summary, paste0(path_output, sample, "_summary.txt"))