library(dplyr)
library(readxl)
library(ggplot2)
library(patchwork)

path_output <- "figures/figure2/"

datasets <- c("LUAD", "mCRPC")
marker <- "CD68"
min_score <- 0.4
cofactor <- 10
colours <- c(neg = "#3973AC", pos = "#D1495B")

# Expected marker-positive and marker-negative cell types (Table S5)
posneg <- read_xlsx("utils/protein_posneg.xlsx") %>% filter(protein == marker)
stopifnot(nrow(posneg) == 1)
positive <- strsplit(posneg$Positive, ";")[[1]]
negative <- strsplit(posneg$Negative, ";")[[1]]

# asinh-normalized marker intensity per cell with its cell type
marker_data <- bind_rows(lapply(datasets, function(dataset) {
  # Raw protein intensities of all cells (columns named <sample>_<cell>)
  protein <- readRDS(paste0("output/protein_expression/prot_", dataset, ".rds"))
  rownames(protein) <- sub("\\.1$", "", rownames(protein))
  stopifnot(marker %in% rownames(protein))
  cell_sample <- sub("_.*", "", colnames(protein))
  
  # Cell type labels, low-confidence and unlabelled cells as Unknown
  labels <- bind_rows(lapply(unique(cell_sample), function(sample) {
    predictions <- read.table(paste0("output/cell_annotation/", dataset, "/", sample, "_celltype.txt"))
    data.frame(cell = paste0(sample, "_", rownames(predictions)),
               cell_type = if_else(predictions$prediction.score.max > min_score, predictions$predicted.id, "Unknown"))
  }))
  cell_type <- labels$cell_type[match(colnames(protein), labels$cell)]
  cell_type[is.na(cell_type)] <- "Unknown"
  
  data.frame(dataset = dataset,
             cell_type = cell_type,
             intensity = asinh(as.numeric(protein[marker, ]) / cofactor))
})) %>%
  mutate(population = case_when(cell_type %in% positive ~ "pos",
                                cell_type %in% negative ~ "neg",
                                TRUE ~ "none"))

zero_stats <- marker_data %>%
  group_by(dataset) %>%
  summarise(pct_zero = mean(intensity == 0) * 100)

# Figure 2D: intensity distribution with expected positive and negative populations
p_main <- ggplot(marker_data, aes(x = intensity)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "grey80", colour = "white") +
  geom_density(data = filter(marker_data, population != "none"),
               aes(colour = population, fill = population), alpha = 0.20, linewidth = 1) +
  geom_text(data = zero_stats, aes(x = 0.25, y = Inf, label = sprintf("Zeros: %.1f%%", pct_zero)),
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.3, size = 3.5) +
  facet_wrap(~dataset, scales = "free", ncol = 1) +
  scale_y_sqrt() +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  labs(x = "Normalized intensity", y = "Density (sqrt scale)", title = marker) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

p_legend <- ggplot() +
  annotate("label", x = 0, y = 1, hjust = 0, vjust = 1, size = 3.4, lineheight = 1.1,
           label = paste0("POSITIVE\n", paste(positive, collapse = "\n")),
           colour = "#B82E45", fill = "#FBE9EC", label.size = 0.3) +
  annotate("label", x = 0, y = 0.48, hjust = 0, vjust = 1, size = 3.4, lineheight = 1.1,
           label = paste0("NEGATIVE\n", paste(negative, collapse = "\n")),
           colour = "#28669C", fill = "#E7F0F8", label.size = 0.3) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(path_output, "fig2d_CD68_posneg.pdf"), p_main + p_legend + plot_layout(widths = c(4, 2)),
       width = 4, height = 5, units = "in")