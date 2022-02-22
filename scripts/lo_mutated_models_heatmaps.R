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
library(rstatix)
library(ggpubr)

# read the parameterization data from the gitsbe link-operator mutated models
if (!file.exists('data/lo_df.rds')) {
  # see Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `parameterization-comp.tar.gz`
  models_dir = '/home/john/tmp/ags-paper/parameterization-comp/link-only/gitsbe_link_only_cascade_2.0_ss_20200805_231150/models'

  lo_df = emba::get_link_operators_from_models_dir(models_dir)

  # save data.frame
  saveRDS(object = lo_df, file = "data/lo_df.rds")
} else {
  lo_df = readRDS(file = "data/lo_df.rds")
}

# remove model names and covert to matrix for the heatmap
lo_mat = as.matrix(lo_df)
rownames(lo_mat) = NULL

# read the stable state data from the gitsbe link-operator mutated models
if (!file.exists('data/lo_ss_df.rds')) {
  # see Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `parameterization-comp.tar.gz`
  models_dir = '/home/john/tmp/ags-paper/parameterization-comp/link-only/gitsbe_link_only_cascade_2.0_ss_20200805_231150/models'

  lo_ss_df = emba::get_stable_state_from_models_dir(models_dir)

  # save data.frame
  saveRDS(object = lo_ss_df, file = "data/lo_ss_df.rds")
} else {
  lo_ss_df = readRDS(file = "data/lo_ss_df.rds")
}

# remove model names and covert to matrix for the heatmap
lo_ss_mat = as.matrix(lo_ss_df)
rownames(lo_ss_mat) = NULL

# get the AGS steady state (TIDY UP THE DATA FROM THE STEADY STATE FILE)
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
  'NF-kB' = '#FFB500', 'Apoptosis' = '#B903AA', 'RTK' = '#B79762',
  'PI3K-AKT' = '#3B5DFF', 'mTOR' = '#00C2A0')

# make a vector of node-pathway annotations
node_path_map = node_path_tbl %>% pull(path_abbrev)
names(node_path_map) = node_path_tbl %>% pull(node)

# get CASCADE 2.0
edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_2.0.sif', delim = " ", col_names = c('source', 'effect', 'target'), col_types = "ccc")

# find number of regulators (target connectivity, in-degree) per node
targets = edge_tbl %>% distinct(target) %>% pull()
node_conn_map = sapply(targets, function(trg) {
  # get number of regulators
  edge_tbl %>% filter(target == trg) %>% distinct(source) %>% summarise(n()) %>% pull()
})

# find source connectivity (out-degree)
sources = edge_tbl %>% distinct(source) %>% pull()
node_src_conn_map = sapply(sources, function(src) {
  edge_tbl %>% filter(source == src) %>% distinct(target) %>% tally() %>% pull()
})

# drug target vector annotation
node_drug_target_map = node_path_tbl %>% pull(is_target)
names(node_drug_target_map) = node_path_tbl %>% pull(node)

# COSMIC annotation
node_cosmic_role = readRDS(file = 'data/cosmic_tbl.rds') # see 'get_cosmic_data_annot.R'
node_cosmic_map = sapply(targets, function(target) { # all 144 CASCADE 2.0 nodes
  role = node_cosmic_role %>%
    filter(cascade2_node == target) %>%
    pull(role)
  if (length(role) > 0 && role == 'oncogene, TSG') role = 'Both' # simplify name
  return(ifelse(length(role) == 0, NA, role))
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

# define coloring for drug target annotation
drug_target_colors = c('FALSE' = 'white', 'TRUE' = 'purple')

# define coloring for the COSMIC annotation
set1_col = RColorBrewer::brewer.pal(9, 'Set1')
cosmic_colors = c('Both' = set1_col[9], 'oncogene' = set1_col[1], 'TSG' = set1_col[3])

# find percent agreement between link-operator parameterization and stable state activity

# check the models are the same (same order in the matrices as well)
stopifnot(all(rownames(lo_df) == rownames(lo_ss_df)))

lo_ss_aggreement = sapply(colnames(lo_mat), function(node) {
  lo_data = lo_mat[,node]
  ss_data = lo_ss_mat[,node]

  add_vec = lo_data + ss_data
  sub_vec = lo_data - ss_data

  # AND-NOT (0) == Inhibited stable state (0)
  and_not_0ss_agreement = sum(add_vec == 0)
  # OR-NOT (1) == Active stable state (1)
  or_not_1ss_agreement = sum(add_vec == 2)
  # AND-NOT (0) != Active stable state (1)
  and_not_1ss_disagreement = sum(sub_vec == -1)
  # OR-NOT  (1) != Inhibited stable state (0)
  or_not_0ss_disagreement = sum(sub_vec == 1)

  n = and_not_0ss_agreement + or_not_1ss_agreement + and_not_1ss_disagreement + or_not_0ss_disagreement
  stopifnot(n == length(ss_data))

  percent_agreement = (and_not_0ss_agreement + or_not_1ss_agreement)/n

  return(percent_agreement)
})

####################################################################
# Compare Oncogenes vs TSGs mean activity values from stable state #
####################################################################

## remove NA's and 'Both' category (oncogenes and TSGs)
node_cosmic_map2 = node_cosmic_map[!is.na(node_cosmic_map) & node_cosmic_map != 'Both']
node_names = node_cosmic_map2 %>% names()
cosmic_annot = node_cosmic_map2 %>% unname()

## Get the average state per node of interest
mean_states = colMeans(lo_ss_mat)[node_names] %>% unname()

cosmic_state = tibble(node = node_names, cosmic = cosmic_annot, mean_state = mean_states)
cosmic_state %>% wilcox_test(mean_state ~ cosmic)

## Apply Wilcox test
stat_test = cosmic_state %>%
  rstatix::wilcox_test(mean_state ~ cosmic) %>%
  rstatix::add_significance("p")

## Visualize data and save plot
ggpubr::ggboxplot(data = cosmic_state %>% rename(COSMIC = cosmic),
  x = 'COSMIC', y = 'mean_state', fill = 'COSMIC', palette = cosmic_colors[2:3],
  xlab = '', ylab = 'Average Activity State', add = "jitter") +
  ylim(c(0,1.15)) +
  ggpubr::stat_pvalue_manual(stat_test, label = "p = {p} ({p.signif})", y.position = 1.1)
ggsave(filename = 'img/cosmic_state_cmp.png', dpi = "print", width = 7, height = 5)
ggsave(filename = 'img/cosmic_state_cmp.pdf', width = 7, height = 5)

#########################
# Link-operator Heatmap #
#########################
# - Column K-means clustering (3)
# - Pathway annotation
# - Connectivity annotation
# - COSMIC annotation

# define annotations
ha_lo = HeatmapAnnotation(COSMIC = node_cosmic_map[colnames(lo_mat)],
  Pathway = node_path_map[colnames(lo_mat)],
  Connectivity = anno_barplot(x = node_conn_map[colnames(lo_mat)]),
  col = list(Pathway = pathway_colors, COSMIC = cosmic_colors),
  annotation_legend_param = list(COSMIC = list(at = c('oncogene', 'TSG', 'Both'))),
  na_col = "white")

# for testing
#indexes = sample(1:nrow(lo_mat), 500)

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
# - COSMIC annotation

# define annotations
ha_ss = HeatmapAnnotation(Training = node_training_state_map, # by default ordered by `colnames(lo_ss_mat)`
  COSMIC = node_cosmic_map[colnames(lo_ss_mat)],
  Pathway = node_path_map[colnames(lo_ss_mat)],
  Connectivity = anno_barplot(x = node_conn_map[colnames(lo_ss_mat)]),
  col = list(Training = training_colors, Pathway = pathway_colors, COSMIC = cosmic_colors),
  na_col = "white",
  annotation_legend_param = list(COSMIC = list(at = c('oncogene', 'TSG', 'Both'))),
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
draw(heatmap_ss, annotation_legend_side = "right", merge_legends = FALSE)
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

##################################################
# Combined Stable States + Link Operator Heatmap #
##################################################
# - Column K-means clustering (3)
# - Training data annotation
# - Pathway annotation
# - Drug Target Annotation
# - Source and Target Connectivity annotations
# - COSMIC annotation
# - Percent Agreement between link-operator and stable state activity

set.seed(42)
heat_ss = Heatmap(matrix = lo_ss_mat[,colnames(lo_mat)],
  #matrix = lo_ss_mat[indexes, colnames(lo_mat)],
  name = 'heat_ss',
  column_km = 3, column_km_repeats = 5, col = state_colors,
  row_title_side = "right", row_title = "Stable states",
  show_row_dend = FALSE, #row_title_rot = 90,
  show_heatmap_legend = TRUE, column_title = "Combined Heatmaps",
  heatmap_legend_param = list(title = 'Activity State', labels = c('Inhibited', 'Active')),
  use_raster = TRUE, raster_quality = 4) # change raster_quality = 1 to avoid having white lines in pdf heatmap

ro_ss = row_order(heat_ss)

heat_param = Heatmap(matrix = lo_mat,
  #matrix = lo_mat[indexes, ],
  name = 'heat_param', column_km = 3, column_km_repeats = 5,
  row_order = ro_ss, # keep same order as in stable state heatmap!
  col = lo_colors, row_title_side = "right", row_title = "Parameterization",
  show_row_dend = FALSE, #row_title_rot = 90,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(title = 'Link Operator', labels = c('AND-NOT', 'OR-NOT')))

ha = HeatmapAnnotation(Calibration = node_training_state_map[colnames(lo_mat)],
  COSMIC = node_cosmic_map[colnames(lo_mat)],
  Pathway = node_path_map[colnames(lo_mat)],
  `Drug Targets` = node_drug_target_map[colnames(lo_mat)],
  `Target Connectivity` = anno_barplot(x = node_conn_map[colnames(lo_mat)]),
  `Source Connectivity` = anno_barplot(x = node_src_conn_map[colnames(lo_mat)]),
  Agreement = anno_barplot(x = lo_ss_aggreement, gp = gpar(fill = 'black')),
  col = list(Calibration = training_colors, Pathway = pathway_colors, COSMIC = cosmic_colors, `Drug Targets` = drug_target_colors),
  na_col = "white",
  annotation_name_side = 'right', annotation_name_rot = list(Connectivity = 0),
  annotation_legend_param = list(COSMIC = list(at = c('oncogene', 'TSG', 'Both'))),
  show_legend = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), gap = unit(c(1,1,1,2,5,5), "points"))

column_name_annot = HeatmapAnnotation(node_names = anno_text(colnames(lo_mat), gp = gpar(fontsize = 6)))

# we combine the heatmaps vertically along with the annotations
heat_list = heat_ss %v% heat_param %v% ha %v% column_name_annot

# mark 3 nodes with boxes
co = column_order(heat_list)$`2` # the marked nodes are in the second column slice
nc = length(co)
marked_nodes = c("JNK_f", "ERK_f", "MAPK14")

png(filename = "img/lo_combined_heat.png", width = 7, height = 7, units = "in", res = 600)
#cairo_pdf(filename = 'img/lo_combined_heat.pdf', width = 7, height = 7) # for pdf output
draw(heat_list, heatmap_legend_side = "left")
decorate_annotation("Agreement", {
  grid.lines(x = unit(c(0.5, 53.5), units = "native"), # 52 link-operator nodes
    y = unit(c(0.5, 0.5), units = "npc"),
    gp = gpar(lty = 2, col = 'red'))
})
# all the marked nodes belong to the second slice
decorate_heatmap_body(heatmap = "heat_ss", column_slice = 2, code = {
  for(node in marked_nodes) {
    i = which(colnames(lo_mat)[co] == node)
    grid.rect(x = (i-0.5)/nc, width = 1/nc, gp=gpar(col="black", fill = NA, lwd = 1))
  }
})
dev.off()

##########################################################
# Reduced Combined Stable States + Link Operator Heatmap #
##########################################################

# Less column annotations will be included (Training and Cosmic)
ha_reduced = HeatmapAnnotation(Calibration = node_training_state_map[colnames(lo_mat)],
  COSMIC = node_cosmic_map[colnames(lo_mat)],
  col = list(Calibration = training_colors, COSMIC = cosmic_colors),
  na_col = "white",
  annotation_name_side = 'right',
  annotation_legend_param = list(COSMIC = list(at = c('oncogene', 'TSG', 'Both'))),
  show_legend = c(FALSE, TRUE), gap = unit(c(1,1), "points"))

# color the 3 marked nodes differently
col_name_colors = rep("black", ncol(lo_mat))
names(col_name_colors) = colnames(lo_mat)
col_name_colors[names(col_name_colors) %in% marked_nodes] = 'blue'

column_name_annot_reduced = HeatmapAnnotation(node_names = anno_text(
  colnames(lo_mat), gp = gpar(fontsize = 8, col = col_name_colors)))

heat_list_reduced = heat_ss %v% heat_param %v% ha_reduced %v% column_name_annot_reduced

png(filename = "img/lo_combined_heat_reduced.png", width = 7, height = 7, units = "in", res = 600)
draw(heat_list_reduced, heatmap_legend_side = "left")
# all the marked nodes belong to the second slice
decorate_heatmap_body(heatmap = "heat_ss", column_slice = 2, code = {
  for(node in marked_nodes) {
    i = which(colnames(lo_mat)[co] == node)
    grid.rect(x = (i-0.5)/nc, width = 1/nc, gp=gpar(col="black", fill = NA, lwd = 2))
  }
})
decorate_heatmap_body(heatmap = "heat_param", column_slice = 2, code = {
  for(node in marked_nodes) {
    i = which(colnames(lo_mat)[co] == node)
    grid.rect(x = (i-0.5)/nc, width = 1/nc, gp=gpar(col="black", fill = NA, lwd = 2))
  }
})
dev.off()
