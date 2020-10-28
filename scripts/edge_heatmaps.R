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

if (!file.exists('data/edge_mat.rds')) {
  # get CASCADE 2.0
  edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_2.0.sif', delim = " ", col_names = c('source', 'effect', 'target'), col_types = "ccc")

  # read gitsbe topology/edge mutated models
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
  saveRDS(object = edge_mat, file = "data/edge_mat.rds")
} else {
  edge_mat = readRDS(file = "data/edge_mat.rds")
}

# heatmap prerequisites

# define coloring for the edges (absence vs presence)
edge_colors = c("red", "green4")
edge_col_fun = circlize::colorRamp2(breaks = c(0,1), colors = edge_colors)

# define the legend
edge_legend = Legend(title = "Edge Mutations",
  labels = c("Absense", "Presence"), legend_gp = gpar(fill = edge_colors))
legend_list = packLegend(edge_legend)

# read node-pathway annotation for CASCADE 2.0
node_path_tbl = readRDS(file = 'data/node_path_tbl.rds')

# Map pathway names to distinct colors
pathway_colors = c('Cross-talk' = 'black', 'MAPK' = 'red',
  'TGF-b' = '#4FC601', 'Wnt' = 'blue', 'Rho' = '#A1C299',
  'Cell Cycle' = '#7A4900', 'JAK-STAT' = '#1CE6FF', 'NF-kB' = '#FF4A46',
  'Apoptosis' = '#B903AA', 'RTK' = '#B79762', 'PI3K-AKT' = '#3B5DFF',
  'mTOR' = '#00C2A0')

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
edge_path_annot = sapply(colnames(edge_mat), function(edge) {
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
stopifnot(all(names(edge_path_annot) == colnames(edge_mat)))

#################################
# Edge Distribution in Pathways #
#################################
dplyr::bind_cols(edge = names(edge_path_annot), pathway = edge_path_annot) %>%
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
pathway_annot = HeatmapAnnotation(Pathway = edge_path_annot,
  name = 'pathway_annot', col = list(Pathway = pathway_colors))

################
# Edge Heatmap #
################
# Extra:
# - Column K-means clustering (4)
# - Pathway Annotation

indx = sample(1:nrow(edge_mat), size = 500)

set.seed(42)
edge_heat = ComplexHeatmap::Heatmap(matrix = edge_mat,
  name = "edge_heatmap", bottom_annotation = pathway_annot,
  column_title = "Model topology parameterization", column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 1), column_km = 4,
  col = edge_col_fun, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = FALSE, use_raster = TRUE,
  raster_device = "png", raster_quality = 20)

png(filename = "img/edge_heat.png", width = 7, height = 5, units = "in", res = 600)
draw(edge_heat, annotation_legend_list = legend_list,
  annotation_legend_side = "right")
dev.off()

###############################################
# Subset Edge Heatmap with K-means clustering #
###############################################
# Extra:
# - Subset columns to the 'stable' edges (those that do not change so much across the models)
# - Column K-means clustering (2)
# - Pathway Annotation

# Subset the matrix data to the 'stable' edges only
edge_avg = edge_mat %>% colSums()/nrow(edge_mat)
stable_edges = names(edge_avg[edge_avg < 0.4 | edge_avg > 0.95]) # user-defined thresholds
stable_edge_mat = edge_mat[,stable_edges]

# Subset the pathway annotation to the 'stable' edges only
stable_edge_path_annot = edge_path_annot[names(edge_path_annot) %in% colnames(stable_edge_mat)]

# sanity data check
stopifnot(all(names(stable_edge_path_annot) == colnames(stable_edge_mat)))

# Define stable edge annotation
pathway_stable_edges_annot = HeatmapAnnotation(Pathway = stable_edge_path_annot,
  name = 'pathway_annot_stable',
  col = list(Pathway = pathway_colors[!names(pathway_colors) %in% c('Cell Cycle', 'PI3K-AKT')]))
# no 'stable' edges in the 'Cell Cycle' and 'PI3K-AKT' pathways!

set.seed(42)
edge_heat_stable = ComplexHeatmap::Heatmap(matrix = stable_edge_mat,
  name = "edge_heatmap", cluster_rows = FALSE,
  bottom_annotation = pathway_stable_edges_annot,
  column_title = "Model topology parameterization (stable edges)",
  column_title_gp = gpar(fontsize = 20),
  column_names_gp = gpar(fontsize = 6), column_km = 2,
  col = edge_col_fun, show_row_names = FALSE, show_row_dend = FALSE,
  show_heatmap_legend = FALSE, use_raster = TRUE, raster_quality = 20)

png(filename = "img/edge_heat_stable.png", width = 7, height = 5, units = "in", res = 600)
draw(edge_heat_stable, annotation_legend_list = legend_list,
  annotation_legend_side = "right")
dev.off()
