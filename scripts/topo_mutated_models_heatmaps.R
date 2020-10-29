#############################################
# Heatmaps for topology/edge mutated models #
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

# read the parameterization data from the gitsbe topology mutated models
if (!file.exists('data/edge_mat.rds')) {
  # get CASCADE 2.0
  edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_2.0.sif', delim = " ", col_names = c('source', 'effect', 'target'), col_types = "ccc")

  # see Zenodo dataset [TOADD], file `parameterization-comp.tar.gz`
  models_dir = '/home/john/tmp/ags-paper/parameterization-comp/topology-only/gitsbe_topology_only_cascade_2.0_ss_20200806_085342/models'

  patterns = edge_tbl %>%
    mutate(pattern_str = paste0(target, ".*", source)) %>%
    pull(pattern_str)

  count = 1
  for (f in list.files(path = models_dir, full.names = TRUE)) {
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
  saveRDS(object = edge_mat, file = "data/edge_mat.rds")
} else {
  edge_mat = readRDS(file = "data/edge_mat.rds")
}

# read the stable state data from the gitsbe topology mutated models
if (!file.exists('data/topo_ss_df.rds')) {
  # see Zenodo dataset [TOADD], file `parameterization-comp.tar.gz`
  models_dir = '/home/john/tmp/ags-paper/parameterization-comp/topology-only/gitsbe_topology_only_cascade_2.0_ss_20200806_085342/models'

  # 70 models has 2 stable states (check with `all.ss = TRUE` parameter)
  topo_ss_df = emba::get_stable_state_from_models_dir(models_dir)

  # save data.frame
  saveRDS(object = topo_ss_df, file = "data/topo_ss_df.rds")
} else {
  topo_ss_df = readRDS(file = "data/topo_ss_df.rds")
}

# remove model names and covert to matrix for the heatmap
rownames(topo_ss_df) = NULL
topo_ss_mat = as.matrix(topo_ss_df)

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

#########################
# heatmap prerequisites #
#########################

# define coloring for the edges (Absence => red vs Presence => green)
edge_colors = c("red", "green4")
names(edge_colors) = c(0,1)

# define coloring for the stable states (Active => red, Inhibited => green)
state_colors = edge_colors

# define the training steady state for node annotation
node_training_state_map = rep(NA, ncol(topo_ss_mat))
names(node_training_state_map) = colnames(topo_ss_mat)

for (node in names(steady_state)) {
  if (unname(steady_state[node]) == 1) {
    node_training_state_map[node] = 'Active'
  } else { # 0
    node_training_state_map[node] = 'Inhibited'
  }
}

training_colors = c('Inhibited' = 'red', 'Active' = 'green4')

# define the legend
# edge_legend = Legend(title = "Edge Mutations",
#   labels = c("Absense", "Presence"), legend_gp = gpar(fill = edge_colors))
# activity_state_legend = Legend(title = "Activity State",
#   labels = c("Inhibited", "Active"), legend_gp = gpar(fill = state_colors))

# read node-pathway annotation for CASCADE 2.0
node_path_tbl = readRDS(file = 'data/node_path_tbl.rds')

# Map pathway names to distinct colors
pathway_colors = c('Cross-talk' = 'black', 'MAPK' = 'red',
  'TGF-b' = '#4FC601', 'Wnt' = 'blue', 'Rho' = '#A1C299',
  'Cell Cycle' = '#7A4900', 'JAK-STAT' = '#1CE6FF', 'NF-kB' = '#FF4A46',
  'Apoptosis' = '#B903AA', 'RTK' = '#B79762', 'PI3K-AKT' = '#3B5DFF',
  'mTOR' = '#00C2A0')

# make a vector of node-pathway annotations
node_path_map = node_path_tbl %>% pull(path_abbrev)
names(node_path_map) = node_path_tbl %>% pull(node)

#################################
# Node Distribution in Pathways #
#################################
node_path_tbl %>%
  group_by(path_abbrev) %>%
  summarize(node_prop = n()/nrow(.), .groups = 'drop') %>%
  rename(pathway = path_abbrev) %>%
  tidyr::drop_na() %>%
  mutate(pathway = forcats::fct_reorder(pathway, desc(node_prop))) %>%
  ggplot(aes(x = pathway, y = node_prop, fill = pathway)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = pathway_colors) +
    scale_y_continuous(labels = scales::percent, limits = c(0,0.3)) +
    labs(x = 'Pathway', y = 'Proportion of total nodes in Pathway',
      title = 'Node Distribution across Pathways in CASCADE 2.0') +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(filename = 'img/node_path_dist.png', dpi = "print", width = 7, height = 5)

# convert to edge-pathway annotation for CASCADE 2.0
edge_path_map = sapply(colnames(edge_mat), function(edge) {
  split_res = stringr::str_split(edge, pattern = ' ', simplify = TRUE)
  source = split_res[1]
  target = split_res[3]
  source_path = node_path_tbl %>% filter(node == source) %>% pull(path_abbrev)
  target_path = node_path_tbl %>% filter(node == target) %>% pull(path_abbrev)

  # in case of logical nodes as edge endpoints (`NA` pathway),
  # return the pathway of the biological node
  if (is.na(source_path)) {
    return(target_path)
  }
  if (is.na(target_path)) {
    return(source_path)
  }

  # return the common pathway name or `Cross-talk` if edge is connecting
  # nodes from different pathways
  if (source_path == target_path) {
    return(source_path)
  } else return('Cross-talk')
})

# sanity data check
stopifnot(all(names(edge_path_map) == colnames(edge_mat)))

#################################
# Edge Distribution in Pathways #
#################################
dplyr::bind_cols(edge = names(edge_path_map), pathway = edge_path_map) %>%
  group_by(pathway) %>%
  summarize(edge_prop = n()/nrow(.), .groups = 'drop') %>%
  mutate(pathway = forcats::fct_reorder(pathway, desc(edge_prop))) %>%
  ggplot(aes(x = pathway, y = edge_prop, fill = pathway)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = pathway_colors) +
    scale_y_continuous(labels = scales::percent, limits = c(0,0.6)) +
    labs(x = 'Pathway', y = 'Proportion of total edges in Pathway',
      title = 'Edge Distribution across Pathways in CASCADE 2.0') +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(filename = 'img/edge_path_dist.png', dpi = "print", width = 7, height = 5)

# define pathway annotation
pathway_annot = HeatmapAnnotation(Pathway = edge_path_map,
  name = 'pathway_annot', col = list(Pathway = pathway_colors))

################
# Edge Heatmap #
################
# Extra:
# - Column K-means clustering (4)
# - Pathway Annotation

indexes = sample(1:nrow(edge_mat), size = 500)

set.seed(42)
edge_heat = ComplexHeatmap::Heatmap(matrix = edge_mat,
  name = "edge_heatmap", bottom_annotation = pathway_annot,
  column_title = "Model topology parameterization", column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 1), column_km = 4,
  col = edge_colors, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Edge Mutations', labels = c('Absense', 'Presence')))
  #,use_raster = TRUE, raster_quality = 20)

#legend_list = packLegend(edge_legend)

png(filename = "img/edge_heat.png", width = 7, height = 5, units = "in", res = 600)
draw(edge_heat, annotation_legend_side = "right", merge_legends = TRUE)
dev.off()

#######################
# Subset Edge Heatmap #
#######################
# Extra:
# - Subset columns to the 'stable' edges (those that do not change so much across the models)
# - Column K-means clustering (2)
# - Pathway Annotation

# Subset the matrix data to the 'stable' edges only
edge_avg = edge_mat %>% colSums()/nrow(edge_mat)
stable_edges = names(edge_avg[edge_avg < 0.4 | edge_avg > 0.95]) # user-defined thresholds
stable_edge_mat = edge_mat[,stable_edges]

# Subset the pathway annotation to the 'stable' edges only
stable_edge_path_map = edge_path_map[names(edge_path_map) %in% colnames(stable_edge_mat)]

# sanity data check
stopifnot(all(names(stable_edge_path_map) == colnames(stable_edge_mat)))

# Define annotations
ha = HeatmapAnnotation(Pathway = stable_edge_path_map,
  col = list(Pathway = pathway_colors[!names(pathway_colors) %in% c('Cell Cycle', 'PI3K-AKT')]))
# no 'stable' edges in the 'Cell Cycle' and 'PI3K-AKT' pathways!

set.seed(42)
edge_heat_stable = ComplexHeatmap::Heatmap(matrix = stable_edge_mat,
  name = "edge_heatmap", cluster_rows = FALSE,
  bottom_annotation = ha,
  column_title = "Model topology parameterization (stable edges)",
  column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 6), column_km = 2,
  col = edge_colors, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Edge Mutations', labels = c('Absense', 'Presence')))
  #, use_raster = TRUE, raster_quality = 20)

#legend_list = packLegend(edge_legend)

png(filename = "img/edge_heat_stable.png", width = 7, height = 5, units = "in", res = 600)
draw(edge_heat_stable, annotation_legend_side = "right", merge_legends = TRUE)
dev.off()

#########################
# Stable States Heatmap #
#########################
# - Column K-means clustering (3)
# - Training data annotation
# - Pathway annotation

# define annotations
node_path_map_ss = node_path_map[colnames(topo_ss_mat)]
# data/order check
stopifnot(all(names(node_path_map_ss) == colnames(topo_ss_mat)))

ha_ss = HeatmapAnnotation(Training = node_training_state_map, Pathway = node_path_map_ss,
  col = list(Training = training_colors, Pathway = pathway_colors),
  na_col = "white",
  show_legend = c("Training" = FALSE))

indexes = sample(1:nrow(topo_ss_mat), size = 500)

set.seed(42)
heatmap_ss = ComplexHeatmap::Heatmap(matrix = topo_ss_mat,
  name = "heatmap_ss", column_km = 3, column_km_repeats = 5,
  bottom_annotation = ha_ss,
  column_title = "Models Stable States", column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 3),
  column_dend_height = unit(1, "inches"),
  col = state_colors, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Activity State', labels = c('Inhibited', 'Active')))
  #, use_raster = TRUE, raster_quality = 20)

#legend_list = packLegend(activity_state_legend)

png(filename = "img/topo_ss_heat.png", width = 7, height = 5, units = "in", res = 600)
draw(heatmap_ss, annotation_legend_side = "right", merge_legends =TRUE)
dev.off()
