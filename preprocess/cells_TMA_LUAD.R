library(arrow)
library(dplyr)
library(sf)

path_raw <- "raw/LUAD/"
path_regions <- " raw/LUAD/TMA_regions/"
path_objects <- "objects/LUAD/"

raw_data <- c("slide1" = "TMA1-SK11-1/",
              "slide2" = "TMA1-SK11-2/")

pixel_size <- 0.2125  # µm per pixel, morphology_focus image

# QuPath core labels (rows O to A, columns 1 to 6) mapped to TMA1-90 on slide1 and TMA91-180 on slide2
core_labels <- paste0(rep(rev(LETTERS[1:15]), each = 6), "-", 1:6)
slide_offset <- c("slide1" = 0, "slide2" = 90)

# Assign each cell centroid to the TMA core polygon containing it
assign_core <- function(cells, tma_regions) {
  points <- st_as_sf(data.frame(cell_id = cells$cell_id,
                                x = cells$x_centroid / pixel_size,
                                y = cells$y_centroid / pixel_size),
                     coords = c("x", "y"), crs = NA)
  joined <- st_join(points, st_set_crs(tma_regions, NA), join = st_within)
  joined <- joined[!duplicated(joined$cell_id), ]
  stopifnot(identical(joined$cell_id, cells$cell_id))
  joined$name
}

cells_tma <- bind_rows(lapply(names(raw_data), function(slide) {
  cells <- read_parquet(paste0(path_raw, raw_data[slide], "cells.parquet"))
  tma_regions <- read_sf(paste0(path_regions, slide, ".geojson"))
  stopifnot(all(tma_regions$name %in% core_labels))

  core <- assign_core(cells, tma_regions)
  core_index <- match(core, core_labels) + slide_offset[slide]
  cells %>% mutate(slide = slide,
                   TMA = if_else(is.na(core_index), NA_character_, paste0("TMA", core_index)))
}))

write_parquet(cells_tma, paste0(path_objects, "cells_TMA.parquet"))
