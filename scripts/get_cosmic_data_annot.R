##############################################################
# CASCADE 2.0 node annotation with COSMIC Cancer Gene Census #
##############################################################
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(httr)
library(jsonlite)
library(ggplot2)
library(forcats)

# read Cancer Gene Census (CGS) COSMIC data from downloaded file
cosmic_file = 'data/cosmic_cancer_gene_census_all_29102020.tsv'
cosmic_data = readr::read_tsv(file = cosmic_file)

# almost all data has Ensemble gene information annotated in `Synonyms` column,
# so we will use that information for cross referencing
# synonyms = cosmic_data %>% pull(Synonyms)
# stringr::str_detect(string = synonyms, pattern = 'ENS')

# read node-pathway annotation for CASCADE 2.0
node_path_tbl = readRDS(file = 'data/node_path_tbl.rds')

data_list = list()
index = 1

cascade2_nodes = node_path_tbl %>% pull(node)
node_num = 1

for (cascade2_node in cascade2_nodes) {
  print(paste0('CASCADE 2.0 node: ', cascade2_node, ' ',
    round(100*node_num/length(cascade2_nodes), digits = 1), '%'))

  hgnc_symbols_str = node_path_tbl %>%
    filter(node == cascade2_node) %>%
    pull(HGNC_symbol)

  # Prosurvival and Antisurvival nodes do not have HGNC symbols
  if (is.na(hgnc_symbols_str)) next

  # might have more than one HGNC symbol (complexes, etc.)
  hgnc_symbols = unlist(stringr::str_split(hgnc_symbols_str, ", "))

  for (hgnc_symbol in hgnc_symbols) {
    # use HGNC REST API
    resp = httr::GET(url = paste0("http://rest.genenames.org/fetch/symbol/", hgnc_symbol),
      config = add_headers(Accept = 'application/json'))

    # see JSON response in a nice format
    # jsonlite::prettify(rawToChar(resp$content))

    # translate response to an R object
    resp_r = jsonlite::fromJSON(rawToChar(resp$content))
    # get the Ensembl ID from the response
    ensembl_gene_id = resp_r$response$docs$ensembl_gene_id

    print(paste0(hgnc_symbol, " - ", ensembl_gene_id))

    # Is the Ensembl ID in the cosmic dataset?
    cosmic_gene_symbol = cosmic_data %>%
      filter(across(Synonyms, ~ grepl(pattern = ensembl_gene_id, x = .))) %>%
      pull(`Gene Symbol`)

    if (length(cosmic_gene_symbol) == 0) next

    # If it is, get the info that we need from the COSMIC dataset
    cosmic_tbl = cosmic_data %>%
      filter(`Gene Symbol` == cosmic_gene_symbol) %>%
      select(COSMIC_gene_symbol = `Gene Symbol`, Tier, Somatic,
        somatic_tumor_type = `Tumour Types(Somatic)`, role = `Role in Cancer`)

    data_list[[index]] = dplyr::bind_cols(cascade2_node = cascade2_node, HGNC_symbols = hgnc_symbols_str, cosmic_tbl)
    index = index + 1
  }
  node_num = node_num + 1
}

res_df = dplyr::bind_rows(data_list)

# Check: 100 matches on 2/11/2020 (otherwise later filtering and
# "hacked" role assignment might not as expected)
#stopifnot(nrow(res_df) == 100)

# Filter data to only stomach/gastric-related tumor types
ags_related_tumor_types = c('diffuse gastric', 'colorectal carcinoma',
  'colorectal', 'small intestine', 'colorectal cancer', 'gastric',
  'colon', 'GIST')
# res_df %>%
#   filter(str_detect(somatic_tumor_type, pattern = paste(ags_related_tumor_types, collapse = "|"))) %>%
#   filter(Tier == 1)

# Remove genes with exclusively 'fusion' role
# Discard the 'fusion' term on the rest of the `role` annotations (keeping only oncogenes, TSG or both)
# Filter data to only the 'Tier 1' class
# Should be 90 rows remaining
res_df_filt = res_df %>%
  filter(role != 'fusion') %>% # genes with only 'fusion' role are discarded
  mutate(role = replace(role, role == 'TSG, fusion', 'TSG')) %>%
  mutate(role = replace(role, role == 'oncogene, fusion', 'oncogene')) %>%
  mutate(role = replace(role, role == 'oncogene, TSG, fusion', 'oncogene, TSG')) %>%
  filter(Tier == 1)

# Get CASCADE 2.0 nodes that have more than 1 match with COSMIC
# data (complexes, families, genes, etc.)
mult_member_nodes = res_df_filt %>%
  group_by(cascade2_node) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(count > 1) %>% pull(cascade2_node)

# Assign a role to these multi-member nodes (mostly using majority rule)
mult_member_nodes_roles = sapply(mult_member_nodes, function(node) {
  roles = res_df_filt %>% filter(cascade2_node == node) %>% pull(role)
  role_tbl = roles %>% table %>% sort(decreasing = TRUE)

  # we only have nodes with members belonging at most to 2 role classes
  if (length(role_tbl) == 1) {
    return(names(role_tbl))
  } else {
    if (role_tbl[1] > role_tbl[2]) {
      return(names(roles_tbl)[1])
    } else { # since it's sorted decreasing, this case is `roles_tbl[1] == roles_tbl[2]`
      return('oncogene, TSG') # kind-of hacked but it works for this dataset
    }
  }
})

# keep only the distinct CASCADE 2.0 nodes
res_df_filt_unique = res_df_filt %>%
  distinct(cascade2_node, .keep_all = TRUE)

# Change the assigned role of the multi-member nodes
for (node in names(mult_member_nodes_roles)) {
  res_df_filt_unique = res_df_filt_unique %>%
    mutate(role = replace(role, cascade2_node == node, mult_member_nodes_roles[node]))
}

# keep only the info that we need (CASCADE 2.0 node name, COSMIC annotated role)
node_cosmic_role = res_df_filt_unique %>% select(cascade2_node, role)

# save object
saveRDS(object = node_cosmic_role, file = 'data/cosmic_tbl.rds')

# Stacked barplot of Role distribution across
node_cosmic_role %>%
  group_by(role) %>%
  summarise(n = n(), prop = n()/nrow(node_cosmic_role), .groups = 'drop') %>%
  mutate(role_col = 'Role') %>%
  mutate(role = forcats::fct_reorder(role, desc(n))) %>%
  ggplot(aes(x = role_col, y = prop, fill = role)) +
    geom_col() +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, 'Set1'),
      name = 'Role in Cancer',
      breaks = c("oncogene", "TSG", "oncogene, TSG"),
      labels = c('Oncogene', 'TSG', 'Both')) +
    geom_text(aes(label = paste0(scales::percent(prop), ' - ', n, ' nodes')),
      vjust = c(-2.3,2.3,-0.5), size = 7) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    labs(x = 'Roles', y = 'Distribution in CASCADE 2.0 nodes',
    title = 'CASCADE 2.0 - COSMIC Cancer Gene Census') +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_blank())
ggsave(filename = 'img/cosmic_cascade2_dist.png', dpi = "print", width = 7, height = 5)

# # Extra
# # Simple filtering of COSMIC data according to the CASCADE 2.0 name (36 matches)
# cancer_gene_data = cosmic_data %>%
#   filter(`Gene Symbol` %in% cascade2_nodes) %>%
#   select(COSMIC_gene_symbol = `Gene Symbol`, Tier, Somatic,
#     somatic_tumor_type = `Tumour Types(Somatic)`, role = `Role in Cancer`)
