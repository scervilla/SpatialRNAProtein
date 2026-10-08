library(Seurat)
library(dplyr)
library(readxl)
library(irr)
library(spdep)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
dataset <- args[2]

path_objects <- paste0("objects/", dataset, "/")
path_output <- paste0("output/concordance/", dataset, "/")

min_counts <- 3
cofactor <- 10
k_neighbors <- 10
n_permutations <- 99
min_pairs <- 3

set.seed(1)

# Protein-gene pairs (proteins without a panel gene are dropped)
pairs <- read_xlsx("utils/protein_gene.xlsx") %>% filter(!is.na(gene))
thresholds <- read.table(paste0("utils/thresholds_protein_", dataset, ".txt"))

object <- readRDS(paste0(path_objects, sample, ".rds"))
object <- subset(object, subset = nCount_Xenium > min_counts)
object <- NormalizeData(object, assay = "Xenium")

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

# Spatial weights: k nearest neighbors on cell centroids
coords <- GetTissueCoordinates(object[["fov"]], which = "centroids")
stopifnot(identical(coords$cell, colnames(object)))
listw <- nb2listw(knn2nb(knearneigh(as.matrix(coords[, c("x", "y")]), k = k_neighbors)),
                  style = "W", zero.policy = TRUE)

# Detection agreement (Cohen's kappa), Spearman correlation and Lee's L per pair
concordance <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  x <- as.numeric(rna[pairs$gene[i], ])
  y <- as.numeric(protein[pairs$protein[i], ])
  rna_pos <- x > 0
  prot_pos <- y > 0
  double_pos <- rna_pos & prot_pos
  lee <- lee.mc(x, y, listw = listw, nsim = n_permutations, zero.policy = TRUE)
  
  data.frame(sample = sample,
             gene = pairs$gene[i],
             protein = pairs$protein[i],
             n_cells = length(x),
             n_rna_pos = sum(rna_pos),
             pct_rna_pos = mean(rna_pos) * 100,
             n_prot_pos = sum(prot_pos),
             pct_prot_pos = mean(prot_pos) * 100,
             n_double_pos = sum(double_pos),
             pct_double_pos = mean(double_pos) * 100,
             detection_kappa = kappa2(data.frame(rna_pos, prot_pos))$value,
             spearman_all = if (length(x) >= min_pairs) cor(x, y, method = "spearman") else NA_real_,
             spearman_pos = if (sum(double_pos) >= min_pairs) cor(x[double_pos], y[double_pos], method = "spearman") else NA_real_,
             lee_L = as.numeric(lee$statistic),
             lee_p = lee$p.value)
}))

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
write.table(concordance, paste0(path_output, sample, "_summary.txt"))