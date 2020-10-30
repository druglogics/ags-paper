#############################################
# Heatmaps for link-operator mutated models #
#############################################

library(dplyr)
library(tibble)
library(readr)
library(stringr)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(ggplot2)
library(forcats)
library(scales)
library(tidyr)

# read the parameterization data from the gitsbe link-operator mutated models
if (!file.exists('data/lo_df.rds')) {
  # see Zenodo dataset [TOADD], file `parameterization-comp.tar.gz`
  models_dir = '/home/john/tmp/ags-paper/parameterization-comp/link-only/gitsbe_link_only_cascade_2.0_ss_20200805_231150/models'

  lo_df = emba::get_link_operators_from_models_dir(models_dir)

  # save data.frame
  saveRDS(object = lo_df, file = "data/lo_df.rds")
} else {
  lo_df = readRDS(file = "data/lo_df.rds")
}

# remove model names and covert to matrix for the heatmap
rownames(lo_df) = NULL
lo_mat = as.matrix(lo_df)

# read the stable state data from the gitsbe link-operator mutated models
if (!file.exists('data/lo_ss_df.rds')) {
  # see Zenodo dataset [TOADD], file `parameterization-comp.tar.gz`
  models_dir = '/home/john/tmp/ags-paper/parameterization-comp/link-only/gitsbe_link_only_cascade_2.0_ss_20200805_231150/models'

  lo_ss_df = emba::get_stable_state_from_models_dir(models_dir)

  # save data.frame
  saveRDS(object = lo_ss_df, file = "data/lo_ss_df.rds")
} else {
  lo_ss_df = readRDS(file = "data/lo_ss_df.rds")
}

# remove model names and covert to matrix for the heatmap
rownames(lo_ss_df) = NULL
lo_ss_mat = as.matrix(lo_ss_df)

# get the AGS steady state
if (!file.exists('data/steady_state.rds')) {
  steady_state_file = 'data/steadystate'
  lines = readLines(steady_state_file)
  ss_data = unlist(strsplit(x = lines[8], split = "\t"))
  ss_mat = stringr::str_split(string = ss_data, pattern = ":", simplify = TRUE)
  colnames(ss_mat) = c("nodes", "states")
  ss_tbl = ss_mat %>% as_tibble() %>% mutate_at(vars(states), as.integer)

  steady_state = ss_tbl %>% pull(states)
  names(steady_state) = ss_tbl %>% pull(nodes)

  saveRDS(object = steady_state, file = 'data/steady_state.rds')
} else {
  steady_state = readRDS(file = 'data/steady_state.rds')
}

# read node-pathway annotation for CASCADE 2.0
node_path_tbl = readRDS(file = 'data/node_path_tbl.rds')

# Map pathway names to distinct colors
pathway_colors = c('MAPK' = 'red', 'TGF-b' = '#4FC601', 'Wnt' = 'blue',
  'Rho' = '#A1C299', 'Cell Cycle' = '#7A4900', 'JAK-STAT' = '#1CE6FF',
  'NF-kB' = '#FF4A46', 'Apoptosis' = '#B903AA', 'RTK' = '#B79762',
  'PI3K-AKT' = '#3B5DFF', 'mTOR' = '#00C2A0')

# make a vector of node-pathway annotations
node_path_map = node_path_tbl %>% pull(path_abbrev)
names(node_path_map) = node_path_tbl %>% pull(node)

# get CASCADE 2.0
edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_2.0.sif', delim = " ", col_names = c('source', 'effect', 'target'), col_types = "ccc")

# find number of regulators (connectivity) per node
targets = edge_tbl %>% distinct(target) %>% pull()
node_conn_map = sapply(targets, function(trg){
  # get number of regulators
  edge_tbl %>% filter(target == trg) %>% distinct(source) %>% summarise(n()) %>% pull()
})

#########################
# heatmap prerequisites #
#########################

# define coloring for the link-operators (AND-NOT => red, OR-NOT => lightyellow)
lo_colors = c("red", "lightyellow")
names(lo_colors) = c(0,1)

# define coloring for the stable states (Active => red, Inhibited => green)
state_colors = c("red", "green4")
names(state_colors) = c(0,1)

# define the training steady state for node annotation
node_training_state_map = rep(NA, ncol(lo_ss_mat))
names(node_training_state_map) = colnames(lo_ss_mat)

for (node in names(steady_state)) {
  if (unname(steady_state[node]) == 1) {
    node_training_state_map[node] = 'Active'
  } else { # 0
    node_training_state_map[node] = 'Inhibited'
  }
}

training_colors = c('Inhibited' = 'red', 'Active' = 'green4')

#########################
# Link-operator Heatmap #
#########################
# - Column K-means clustering (3)
# - Pathway annotation
# - Connectivity annotation

# define annotations
node_path_map_lo = node_path_map[colnames(lo_mat)]

ha_lo = HeatmapAnnotation(Pathway = node_path_map_lo,
  Connectivity = anno_barplot(x = node_conn_map[colnames(lo_mat)]),
  col = list(Pathway = pathway_colors),
  na_col = "black")

indexes = sample(1:nrow(lo_mat), 500)

set.seed(42)
heatmap_param = ComplexHeatmap::Heatmap(matrix = lo_mat,
  #matrix = lo_mat[indexes, ], # take a subset for testing
  name = "heatmap_param", column_km = 3, column_km_repeats = 5,
  bottom_annotation = ha_lo,
  #clustering_distance_rows = 'binary', clustering_distance_columns = 'binary',
  column_title = "Model Link-Operator Parameterization", column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 8),
  col = lo_colors, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  column_dend_height = unit(0.5, "inches"),
  heatmap_legend_param = list(title = 'Link Operator', labels = c('AND-NOT', 'OR-NOT')))
  #, use_raster = TRUE, raster_quality = 20)

png(filename = "img/lo_heat.png", width = 7, height = 5, units = "in", res = 600)
draw(heatmap_param, annotation_legend_side = "right", merge_legends = FALSE)
dev.off()

#########################
# Stable States Heatmap #
#########################
# - Column K-means clustering (3)
# - Training data annotation
# - Pathway annotation
# - Connectivity annotation

# define annotations
node_path_map_ss = node_path_map[colnames(lo_ss_mat)]

ha_ss = HeatmapAnnotation(Training = node_training_state_map, Pathway = node_path_map_ss,
  Connectivity = anno_barplot(x = node_conn_map[colnames(lo_ss_mat)]),
  col = list(Training = training_colors, Pathway = pathway_colors),
  na_col = "white",
  show_legend = c("Training" = FALSE))

set.seed(42)
heatmap_ss = ComplexHeatmap::Heatmap(matrix = lo_ss_mat,
  #matrix = lo_ss_mat[indexes, ], # take a subset for testing
  name = "heatmap_ss", column_km = 3, column_km_repeats = 5,
  bottom_annotation = ha_ss,
  column_title = "Models Stable States", column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 3),
  column_dend_height = unit(1, "inches"),
  col = state_colors, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Activity State', labels = c('Inhibited', 'Active')))
  #,use_raster = TRUE, raster_quality = 20)

png(filename = "img/lo_ss_heat.png", width = 7, height = 5, units = "in", res = 600)
draw(heatmap_ss, annotation_legend_side = "right", merge_legends = TRUE)
dev.off()

# code to add rectangular boxes for specific nodes

# co = column_order(heatmap_ss)
# nc = ncol(lo_ss_mat)
# marked_nodes = c("MYC", "TP53")
# decorate_heatmap_body(heatmap = "heatmap_ss", code = {
#   for(node in marked_nodes) {
#     i = which(colnames(lo_ss_mat)[co] == node)
#     grid.rect(x = (i-0.5)/nc, width = 1/nc, gp=gpar(col="black", fill = NA, lwd = 1))
#   }
# })
