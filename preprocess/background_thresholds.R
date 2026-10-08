library(Seurat)
library(dplyr)
library(tidyr)
library(readxl)

args = commandArgs(trailingOnly=TRUE)
dataset <- args[1]

path_objects <- c("LUAD" = "objects/LUAD/", "mCRPC" = "objects/mCRPC/")[[dataset]]
path_annotation <- paste0("output/cell_annotation/", dataset, "/")
path_output <- paste0("utils/")

min_counts <- 3
min_score <- 0.4
min_cells <- 5
  
samples <- readLines(paste0("utils/samples_", dataset, ".txt"))

# Marker-negative cell types per protein (Table S5)
negative <- read_xlsx("utils/protein_posneg.xlsx") %>%
  select(protein, cell_type = Negative) %>%
  separate_rows(cell_type, sep = ";")

# Per-sample threshold: 75th percentile of the 90th percentile of each marker-negative cell type
thresholds_sample <- bind_rows(lapply(samples, function(sample) {
  object <- readRDS(paste0(path_objects, sample, ".rds"))
  object <- subset(object, subset = nCount_Xenium > min_counts)
  predictions <- read.table(paste0(path_annotation, sample, "_celltype.txt"))
  stopifnot(identical(rownames(predictions), colnames(object)))
  
  # Protein features share names with genes and are suffixed ".1" by Seurat
  protein <- object[["Protein"]]$counts
  rownames(protein) <- sub("\\.1$", "", rownames(protein))
  stopifnot(all(negative$protein %in% rownames(protein)))
  
  as.data.frame(t(as.matrix(protein[unique(negative$protein), ]))) %>%
    mutate(cell_type = if_else(predictions$prediction.score.max < min_score, "Unknown", predictions$predicted.id)) %>%
    pivot_longer(-cell_type, names_to = "protein", values_to = "intensity") %>%
    inner_join(negative, by = c("protein", "cell_type")) %>%
    group_by(protein, cell_type) %>%
    filter(n() > min_cells) %>%
    summarise(q90 = quantile(intensity, 0.9, na.rm = TRUE), .groups = "drop") %>%
    group_by(protein) %>%
    summarise(threshold = quantile(q90, 0.75)) %>%
    complete(protein = unique(negative$protein)) %>%
    mutate(sample = sample)
}))

# Final threshold: median across samples
thresholds <- thresholds_sample %>%
  group_by(protein) %>%
  summarise(threshold = median(threshold))
stopifnot(!anyNA(thresholds$threshold))

write.table(thresholds_sample, paste0(path_output, "thresholds_protein_", dataset, "_sample.txt"))
write.table(thresholds, paste0(path_output, "thresholds_protein_", dataset, ".txt"))