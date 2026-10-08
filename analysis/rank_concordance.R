library(dplyr)
library(tidyr)
library(readxl)

path_output <- "output/rank_concordance/"

datasets <- c("LUAD", "mCRPC")
metrics <- c("mean_expression", "W0", "W1", "W2", "beta")
spatial_metrics <- c("W0", "W1", "W2", "beta")
min_samples <- 3

# Protein-gene pairs (proteins without a panel gene are dropped)
pairs <- read_xlsx("utils/protein_gene.xlsx") %>% filter(!is.na(gene))

# Per-sample Minkowski profiles (from minkowski_metrics.R)
profiles <- bind_rows(lapply(datasets, function(dataset) {
  samples <- readLines(paste0("utils/samples_", dataset, ".txt"))
  bind_rows(lapply(samples, function(sample) {
    read.table(paste0("output/minkowski/", dataset, "/", sample, ".txt"))
  })) %>%
    mutate(dataset = dataset)
}))

# Paired RNA and protein values per sample, pair and metric;
# a gene backing several proteins (PTPRC: CD45, CD45RA, CD45RO) is copied to each pair
paired <- bind_rows(
  profiles %>%
    filter(assay == "Protein") %>%
    inner_join(pairs, by = c("feature" = "protein")) %>%
    mutate(pair = feature),
  profiles %>%
    filter(assay == "Xenium") %>%
    inner_join(pairs, by = c("feature" = "gene"), relationship = "many-to-many") %>%
    mutate(pair = protein)
) %>%
  select(sample, dataset, pair, assay, all_of(metrics)) %>%
  pivot_longer(all_of(metrics), names_to = "metric") %>%
  pivot_wider(names_from = assay, values_from = value) %>%
  rename(rna = Xenium, protein = Protein)

# Kendall's tau-b between RNA and protein across samples, per dataset, pair and metric
kendall <- paired %>%
  group_by(dataset, pair, metric) %>%
  summarise(n = sum(complete.cases(rna, protein)),
            kendall = if (n >= min_samples) cor(rna, protein, method = "kendall", use = "complete.obs") else NA_real_,
            .groups = "drop")

# Paired Wilcoxon signed-rank test of each spatial metric against mean expression across pairs (both datasets)
kendall_wide <- kendall %>%
  select(dataset, pair, metric, kendall) %>%
  pivot_wider(names_from = metric, values_from = kendall)
wilcoxon <- bind_rows(lapply(spatial_metrics, function(metric) {
  test <- wilcox.test(kendall_wide$mean_expression, kendall_wide[[metric]], paired = TRUE)
  data.frame(comparison = paste0("mean_expression vs ", metric),
             statistic = unname(test$statistic),
             p_value = test$p.value)
})) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
write.table(paired, paste0(path_output, "minkowski_paired.txt"))
write.table(kendall, paste0(path_output, "kendall.txt"))
write.table(wilcoxon, paste0(path_output, "wilcoxon.txt"))
