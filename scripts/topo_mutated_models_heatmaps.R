#############################################
# Heatmaps for topology/edge mutated models #
#############################################

library(dplyr)
library(tibble)
library(readr)
library(stringr)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(ggplot2)
library(forcats)
library(scales)
library(tidyr)

# read the parameterization data from the gitsbe topology mutated models
if (!file.exists('data/edge_mat.rds')) {
  # get CASCADE 2.0
  edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_2.0.sif', delim = " ", col_names = c('source', 'effect', 'target'), col_types = "ccc")

  # see Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `parameterization-comp.tar.gz`
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
  # see Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `parameterization-comp.tar.gz`
  models_dir = '/home/john/tmp/ags-paper/parameterization-comp/topology-only/gitsbe_topology_only_cascade_2.0_ss_20200806_085342/models'

  # 70 models has 2 stable states (check with `all.ss = TRUE` parameter)
  topo_ss_df = emba::get_stable_state_from_models_dir(models_dir)

  # save data.frame
  saveRDS(object = topo_ss_df, file = "data/topo_ss_df.rds")
} else {
  topo_ss_df = readRDS(file = "data/topo_ss_df.rds")
}

# remove model names and covert to matrix for the heatmap
topo_ss_mat = as.matrix(topo_ss_df)
rownames(topo_ss_mat) = NULL

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

# read node-pathway annotation for CASCADE 2.0
node_path_tbl = readRDS(file = 'data/node_path_tbl.rds')

# get CASCADE 2.0
edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_2.0.sif', delim = " ", col_names = c('source', 'effect', 'target'), col_types = "ccc")

# find number of regulators (connectivity, in-degree) per node
targets = edge_tbl %>% distinct(target) %>% pull()
node_conn_map = sapply(targets, function(trg){
  # get number of regulators
  edge_tbl %>% filter(target == trg) %>% distinct(source) %>% summarise(n()) %>% pull()
})

# find source connectivity (out-degree) per node
sources = edge_tbl %>% distinct(source) %>% pull()
node_src_conn_map = sapply(sources, function(src) {
  edge_tbl %>% filter(source == src) %>% distinct(target) %>% tally() %>% pull()
})

# make a vector of edge-target in-degree connectivity annotations
edge_conn_map = sapply(colnames(edge_mat), function(edge) {
  split_res = stringr::str_split(edge, pattern = ' ', simplify = TRUE)
  target = split_res[3]
  unname(node_conn_map[target]) # in-degree of edge's target
})

# make a vector of edge-source out-degree connectivity annotations
edge_src_conn_map = sapply(colnames(edge_mat), function(edge) {
  split_res = stringr::str_split(edge, pattern = ' ', simplify = TRUE)
  source = split_res[1]
  unname(node_src_conn_map[source]) # out-degree or edge's source
})

# Map pathway names to distinct colors
pathway_colors = c('Cross-talk' = 'black', 'MAPK' = 'red',
  'TGF-b' = '#4FC601', 'Wnt' = 'blue', 'Rho' = '#A1C299',
  'Cell Cycle' = '#7A4900', 'JAK-STAT' = '#1CE6FF', 'NF-kB' = '#FFB500',
  'Apoptosis' = '#B903AA', 'RTK' = '#B79762', 'PI3K-AKT' = '#3B5DFF',
  'mTOR' = '#00C2A0')

# make a vector of node-pathway annotations
node_path_map = node_path_tbl %>% pull(var = path_abbrev, name = node)

# drug target edge annotation (drug target can be the source of an edge,
# the target, both source and target or none of them!)
node_drug_target_map = node_path_tbl %>% pull(is_target)
names(node_drug_target_map) = node_path_tbl %>% pull(node)

edge_drug_target_map = sapply(colnames(edge_mat), function(edge) {
  split_res = stringr::str_split(edge, pattern = ' ', simplify = TRUE)
  source = split_res[1]
  target = split_res[3]
  is_edge_src_drug_target = unname(node_drug_target_map[source])
  is_edge_trg_drug_target = unname(node_drug_target_map[target])
  if (is_edge_src_drug_target & is_edge_trg_drug_target)
    return('both')
  if (is_edge_src_drug_target & !is_edge_trg_drug_target)
    return('source')
  if (!is_edge_src_drug_target & is_edge_trg_drug_target)
    return('target')
  if (!is_edge_src_drug_target & !is_edge_trg_drug_target)
    return('none')
})

# colors for edge-drug target annotation
set1_col = RColorBrewer::brewer.pal(9, 'Set1')
drug_target_colors = c('both' = set1_col[2], 'source' = set1_col[1],
  'target' = set1_col[3], 'none' = 'black')

# COSMIC annotation
node_cosmic_role = readRDS(file = 'data/cosmic_tbl.rds') # see 'get_cosmic_data_annot.R'
node_cosmic_map = sapply(targets, function(target) { # all 144 CASCADE 2.0 nodes
  role = node_cosmic_role %>%
    filter(cascade2_node == target) %>%
    pull(role)
  if (length(role) > 0 && role == 'oncogene, TSG') role = 'Both' # simplify name
  return(ifelse(length(role) == 0, NA, role))
})

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
  # or return an 'undirected' merged version of the two pathway names
  # } else return(paste(sort(c(source_path,target_path)), collapse = ","))
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

################
# Edge Heatmap #
################
# Extra:
# - Column K-means clustering (4)
# - Pathway Annotation
# - Drug Target Annotation
# - Target Node In-degree Connectivity
# - Source Node Out-degree Connectivity

# define annotations
ha_edges = HeatmapAnnotation(Pathway = edge_path_map[colnames(edge_mat)],
  `Drug Target` = edge_drug_target_map[colnames(edge_mat)],
  `Target In-degree` = anno_barplot(x = edge_conn_map[colnames(edge_mat)]),
  `Source Out-degree` = anno_barplot(x = edge_src_conn_map[colnames(edge_mat)]),
  col = list(Pathway = pathway_colors, `Drug Target` = drug_target_colors),
  gap = unit(c(1,1,5), "points"))

indexes = sample(1:nrow(edge_mat), size = 500)

set.seed(42)
cl = kmeans(t(edge_mat), centers = 4)$cluster

# hack: split the 3rd cluster to 2 separate ones, in order to distinguish
# the edges whose target has in-degree 1 (a single regulator)
tmp_edge_trg_reg = edge_conn_map[names(cl)[cl == 3]]
edges_with_single_reg_targets = names(tmp_edge_trg_reg)[tmp_edge_trg_reg == 1]
cl[edges_with_single_reg_targets] = 5

edge_heat = ComplexHeatmap::Heatmap(matrix = edge_mat,
  #matrix = edge_mat[indexes, ], # take a subset for testing
  name = "edge_heatmap", bottom_annotation = ha_edges,
  column_title = "Model topology parameterization", column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 1),
  column_split = factor(cl, levels = c(5,3,2,1,4)), cluster_column_slices = FALSE,
  col = edge_colors, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Edge Mutations', labels = c('Absense', 'Presence')),
  use_raster = TRUE, raster_quality = 4)

png(filename = "img/edge_heat.png", width = 7, height = 5, units = "in", res = 600)
draw(edge_heat, annotation_legend_side = "right", merge_legends = FALSE)
dev.off()

#######################
# Subset Edge Heatmap #
#######################
# Extra:
# - Subset columns to the 'stable' edges (those that do not change so much across the models)
# - Column K-means clustering (2)
# - Pathway Annotation
# - Target Node Connectivity
# - Drug Target Annotation
# - Target Node In-degree Connectivity
# - Source Node Out-degree Connectivity

# Subset the matrix data to the 'stable' edges only
edge_avg = edge_mat %>% colSums()/nrow(edge_mat)
# 0.99 => change to 1 if you want to keep all the edges that are always there (low connectivity)
stable_edges = names(edge_avg[edge_avg < 0.45 | (edge_avg > 0.87 & edge_avg <= 0.99)]) # user-defined thresholds
stable_edge_mat = edge_mat[,stable_edges]

# Subset the pathway annotation to the 'stable' edges only
stable_edge_path_map = edge_path_map[names(edge_path_map) %in% colnames(stable_edge_mat)]

# sanity data check
stopifnot(all(names(stable_edge_path_map) == colnames(stable_edge_mat)))

# Define annotations
ha_edges_stable = HeatmapAnnotation(Pathway = stable_edge_path_map,
  `Drug Target` = edge_drug_target_map[colnames(stable_edge_mat)],
  `Target In-degree` = anno_barplot(x = edge_conn_map[colnames(stable_edge_mat)]),
  `Source Out-degree` = anno_barplot(x = edge_src_conn_map[colnames(stable_edge_mat)]),
  annotation_name_gp = gpar(fontsize = c(12,12,11,12)),
  col = list(Pathway = pathway_colors, `Drug Target` = drug_target_colors),
  gap = unit(c(1,2,5), "points"))

set.seed(42)
edge_heat_stable = ComplexHeatmap::Heatmap(matrix = stable_edge_mat,
  #matrix = stable_edge_mat[indexes, ], # take a subset for testing
  name = "edge_heatmap", cluster_rows = FALSE,
  bottom_annotation = ha_edges_stable,
  column_title = "Model topology parameterization (stable edges)",
  column_title_gp = gpar(fontsize = 17),
  column_names_gp = gpar(fontsize = 6), column_km = 2,
  col = edge_colors, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Edge Mutations', labels = c('Absense', 'Presence')))
  #use_raster = TRUE, raster_quality = 20)

png(filename = "img/edge_heat_stable.png", width = 7, height = 5, units = "in", res = 600)
draw(edge_heat_stable, annotation_legend_side = "right", merge_legends = FALSE)
dev.off()

#########################
# Stable States Heatmap #
#########################
# - Column K-means clustering (3)
# - Training data annotation
# - Pathway annotation
# - Connectivity annotation
# - COSMIC annotation

# define coloring for the COSMIC annotation
cosmic_colors = c('Both' = set1_col[4], 'oncogene' = set1_col[1], 'TSG' = set1_col[7])

# define annotations
ha_ss = HeatmapAnnotation(Training = node_training_state_map[colnames(topo_ss_mat)],
  COSMIC = node_cosmic_map[colnames(topo_ss_mat)],
  Pathway = node_path_map[colnames(topo_ss_mat)],
  Connectivity = anno_barplot(x = node_conn_map[colnames(topo_ss_mat)]),
  col = list(Training = training_colors, Pathway = pathway_colors, COSMIC = cosmic_colors),
  na_col = "white",
  annotation_legend_param = list(COSMIC = list(at = c('TSG', 'oncogene', 'Both'))),
  show_legend = c("Training" = FALSE))

indexes = sample(1:nrow(topo_ss_mat), size = 500)

set.seed(42)
heatmap_ss = ComplexHeatmap::Heatmap(matrix = topo_ss_mat,
  #matrix = topo_ss_mat[indexes, ], # take a subset for testing
  name = "heatmap_ss", column_km = 3, column_km_repeats = 5,
  bottom_annotation = ha_ss,
  column_title = "Models Stable States", column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 3),
  column_dend_height = unit(1, "inches"),
  col = state_colors, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Activity State', labels = c('Inhibited', 'Active')))
  #, use_raster = TRUE, raster_quality = 20)

png(filename = "img/topo_ss_heat.png", width = 7, height = 5, units = "in", res = 600)
draw(heatmap_ss, annotation_legend_side = "right", merge_legends = FALSE)
dev.off()
