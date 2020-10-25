library(dplyr)
library(tibble)
library(readr)
library(stringr)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

if (!file.exists('results/edge_mat.rds')) {
  # get CASCADE 2.0
  edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_2.0.sif', delim = " ", col_names = c('source', 'effect', 'target'), col_types = "ccc")

  # read gitsbe topology mutated models
  # see Zenodo dataset [TOADD], file `parameterization-comp.tar.gz`
  data_dir = '/home/john/tmp/ags-paper/parameterization-comp/topology-only/gitsbe_topology_only_cascade_2.0_ss_20200806_085342/models'

  patterns = edge_tbl %>%
    mutate(pattern_str = paste0(target, ".*", source)) %>%
    pull(pattern_str)

  count = 1
  for (f in list.files(path = data_dir, full.names = TRUE)) {
    model_name = stringr::str_remove(string = basename(f), pattern = ".gitsbe")
    print(paste0(model_name, " - ", count))
    count = count + 1

    lines = readr::read_lines(file = f)
    lines = stringr::str_subset(string = lines, pattern = "equation:")

    edges = sapply(patterns, function(edge_pattern) {
      matches = stringr::str_detect(string = lines, pattern = edge_pattern)
      ifelse(test = sum(matches) == 1, TRUE, FALSE)
    }, USE.NAMES = FALSE)

    edge_tbl = edge_tbl %>% tibble::add_column(!!as.symbol(model_name) := edges %>% as.integer())
  }

  # make an edge matrix
  # values denote presence or absence of an edge in the boolean models
  edge_annot = edge_tbl %>%
    mutate(edge = paste(source, effect, target)) %>% pull(edge)

  edge_mat = edge_tbl %>%
    select(-all_of(c("source", "effect", "target"))) %>% t()

  rownames(edge_mat) = NULL
  colnames(edge_mat) = edge_annot

  # save matrix
  saveRDS(object = edge_mat, file = "results/edge_mat.rds")
} else {
  edge_mat = readRDS(file = "results/edge_mat.rds")
}

# heatmap prerequisites

# define coloring
edge_colors = c("red", "green4")
edge_col_fun = circlize::colorRamp2(breaks = c(0, 1), colors = edge_colors)

# define the legends
edge_legend = Legend(title = "Edge Mutations", nrow = 1,
  labels = c("Absense", "Presence"), legend_gp = gpar(fill = edge_colors))
legend_list = packLegend(edge_legend, direction = "vertical")

################
# Edge Heatmap #
################

set.seed(42)
edge_heat = ComplexHeatmap::Heatmap(matrix = edge_mat,
  name = "edge_heatmap",
  column_title = "Model topology parameterization", column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 1),
  col = edge_col_fun, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = FALSE, use_raster = TRUE,
  raster_device = "png", raster_quality = 20)

png(filename = "img/edge_heat.png", width = 7, height = 5, units = "in", res = 600)
draw(edge_heat, annotation_legend_list = legend_list,
  annotation_legend_side = "bottom")
dev.off()

############################################
# Edge Heatmap (column K-means clustering) #
############################################

set.seed(42)
edge_heat = ComplexHeatmap::Heatmap(matrix = edge_mat,
  name = "edge_heatmap",
  column_title = "Model topology parameterization", column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 1), column_km = 4,
  col = edge_col_fun, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = FALSE, use_raster = TRUE,
  raster_device = "png", raster_quality = 20)

png(filename = "img/edge_heat_km.png", width = 7, height = 5, units = "in", res = 600)
draw(edge_heat, annotation_legend_list = legend_list,
  annotation_legend_side = "bottom")
dev.off()

###############################################
# Subset Edge Heatmap with K-means clustering #
###############################################

# keep edges if they do not change too much across the models ('stable' parameterization)
edge_avg = edge_mat %>% colSums()/nrow(edge_mat)

# Subset the matrix to the 'stable' edges only
stable_edges = names(edge_avg[edge_avg < 0.3 | edge_avg > 0.95])
stable_edge_mat = edge_mat[,stable_edges]

set.seed(42)
edge_heat_stable = ComplexHeatmap::Heatmap(matrix = stable_edge_mat,
  name = "edge_heatmap", cluster_rows = FALSE,
  column_title = "Model topology parameterization (stable edges)",
  column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 7), column_km = 2,
  col = edge_col_fun, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = FALSE, use_raster = TRUE, raster_quality = 20)

png(filename = "img/edge_heat_stable_km.png", width = 7, height = 5, units = "in", res = 600)
draw(edge_heat_stable, annotation_legend_list = legend_list,
  annotation_legend_side = "bottom")
dev.off()
