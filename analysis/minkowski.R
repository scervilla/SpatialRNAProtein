library(Seurat)
library(dplyr)
library(readxl)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
dataset <- args[2]

path_objects <- paste0("objects/", dataset, "/")
path_output <- paste0("output/minkowski/", dataset, "/")

min_counts <- 3
cofactor <- 10
bin_size <- 20  # µm
level <- 1      # binarization threshold on mean expression per bin

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

# Features x cells, ordered as the centroids
coords <- GetTissueCoordinates(object[["fov"]], which = "centroids")
genes <- unique(pairs$gene)
expression <- rbind(protein[pairs$protein, coords$cell],
                    as.matrix(rna[genes, coords$cell]))
assay <- rep(c("Protein", "Xenium"), c(nrow(pairs), length(genes)))

# Square grid per core
x_min <- floor(min(coords$x) / bin_size) * bin_size
y_min <- floor(min(coords$y) / bin_size) * bin_size
x_edges <- x_min + 0:max(1, ceiling((max(coords$x) - x_min) / bin_size)) * bin_size
y_edges <- y_min + 0:max(1, ceiling((max(coords$y) - y_min) / bin_size)) * bin_size
n_x <- length(x_edges) - 1
n_y <- length(y_edges) - 1
x_centres <- head(x_edges, -1) + bin_size / 2
y_centres <- head(y_edges, -1) + bin_size / 2

ix <- cut(coords$x, x_edges, include.lowest = TRUE, labels = FALSE)
iy <- cut(coords$y, y_edges, include.lowest = TRUE, labels = FALSE)
stopifnot(!anyNA(ix), !anyNA(iy))
bin <- ix + (iy - 1) * n_x

# Tissue mask: bins with at least one cell
cells_per_bin <- tabulate(bin, nbins = n_x * n_y)
mask <- matrix(cells_per_bin >= 1, n_x, n_y)
n_tissue <- sum(mask)

# Mean expression per bin (empty bins = 0)
membership <- sparseMatrix(i = seq_along(bin), j = bin, x = 1, dims = c(length(bin), n_x * n_y))
bin_mean <- sweep(as.matrix(expression %*% membership), 2, pmax(cells_per_bin, 1), "/")

# 2x2 pixel configurations; minus-sampling keeps blocks with all four corners in tissue
inner_r <- 2:(n_x + 1); inner_c <- 2:(n_y + 1)
r1 <- 1:(n_x + 1); r2 <- 2:(n_x + 2)
c1 <- 1:(n_y + 1); c2 <- 2:(n_y + 2)
mask_p <- matrix(FALSE, n_x + 2, n_y + 2)
mask_p[inner_r, inner_c] <- mask
valid <- mask_p[r1, c1] & mask_p[r1, c2] & mask_p[r2, c1] & mask_p[r2, c2]

# Lookup by number of foreground corners (Michielsen & De Raedt 2001)
area_lut <- c(0, 0.25, 0.5, 0.75, 1)
perim_lut <- c(0, 1, 1, 1, 0)
chi_lut <- c(0, 0.25, 0, -0.25, 0)

results <- bind_rows(lapply(seq_len(nrow(expression)), function(f) {
  density <- matrix(bin_mean[f, ], n_x, n_y)
  
  # Scalar Minkowski functionals W0 (area), W1 (perimeter), W2 (Euler characteristic)
  field_p <- matrix(0, n_x + 2, n_y + 2)
  field_p[inner_r, inner_c] <- density >= level
  TL <- field_p[r1, c1]; TR <- field_p[r1, c2]
  BL <- field_p[r2, c1]; BR <- field_p[r2, c2]
  count <- TL + TR + BL + BR
  diagonal <- (TL == BR) & (TR == BL) & (TL != TR)
  area <- area_lut[count + 1]
  perim <- perim_lut[count + 1]; perim[diagonal] <- 2
  chi <- chi_lut[count + 1]; chi[diagonal] <- -0.5
  
  # Minkowski boundary tensor W1^{0,2} along the isocontour, non-tissue bins masked
  z <- density
  z[!mask] <- NA
  contours <- contourLines(x_centres, y_centres, z, levels = level)
  dx <- unlist(lapply(contours, function(co) diff(co$x)))
  dy <- unlist(lapply(contours, function(co) diff(co$y)))
  ds <- sqrt(dx^2 + dy^2)
  keep <- is.finite(ds) & ds > 0
  dx <- dx[keep]; dy <- dy[keep]; ds <- ds[keep]
  nx <- -dy / ds; ny <- dx / ds  # unit normal = tangent rotated 90 degrees
  
  tensor <- 0.5 * matrix(c(sum(nx * nx * ds), sum(nx * ny * ds),
                           sum(nx * ny * ds), sum(ny * ny * ds)), 2, 2)
  eig <- eigen(tensor, symmetric = TRUE)
  normal_angle <- (atan2(eig$vectors[2, 1], eig$vectors[1, 1]) * 180 / pi) %% 180
  tensor_stats <- c(Txx = tensor[1, 1], Txy = tensor[1, 2], Tyy = tensor[2, 2],
                    lambda_max = eig$values[1], lambda_min = eig$values[2],
                    beta = if (eig$values[1] > 0) eig$values[2] / eig$values[1] else NA_real_,
                    normal_angle = normal_angle,
                    boundary_angle = (normal_angle + 90) %% 180)
  if (length(ds) == 0) tensor_stats[] <- NA_real_
  
  data.frame(level = level,
             W0 = sum(area[valid]) / n_tissue,
             W1 = sum(perim[valid]) / n_tissue,
             W2 = sum(chi[valid]) / n_tissue,
             as.list(tensor_stats),
             contour_length = sum(ds),
             feature = rownames(expression)[f],
             assay = assay[f],
             mean_expression = mean(expression[f, ]),
             pct_expressing = 100 * mean(expression[f, ] > 0),
             sample = sample)
}))

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
write.table(results, paste0(path_output, sample, ".txt"))