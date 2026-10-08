library(Seurat)
library(ggplot2)
library(scCustomize)

path_output <- "figures/figure2/"

patient <- 1
features <- c("Ki-67", "MKI67", "STMN1", "CDK1")

# Patient object with UMAP and normalized RNA and protein (from proliferative_cluster.R)
object <- readRDS(paste0("objects/LUAD/patient", patient, ".rds"))
stopifnot(all(features %in% c(rownames(object[["Xenium"]]), rownames(object[["Protein"]]))))

# Figure 2H: Ki-67 protein and proliferation genes on the UMAP
p <- FeaturePlot_scCustom(object, features = features, order = TRUE, pt.size = 0.001, figure_plot = TRUE)

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig2h_proliferation_umap.pdf"), p, width = 7, height = 7, units = "in")