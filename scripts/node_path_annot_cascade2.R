library(dplyr)
library(readr)
library(stringr)

# read CASCADE 2.0 pathway annotation file
node_path_tbl = readr::read_delim(file = 'data/node_pathway_annotations_cascade2.csv', delim = ";")

# same content columns
stopifnot(all(node_path_tbl$canonicalName == node_path_tbl$name))
stopifnot(all(node_path_tbl$canonicalName == node_path_tbl$`shared name`))

# drop unnecessary/redundant columns
node_path_tbl = node_path_tbl %>%
  select(-one_of(c("canonicalName", "SUID", "hiddenLabel", "selected", "shared name"))) %>%
  relocate(name)

# rename the rest for clarity
node_path_tbl = node_path_tbl %>%
  rename(node = name, is_target = Drug.target,
    pathway = Generic.Pathway.manual,
    HGNC_symbol = HGNC.symbols, in_cascade1 = in_AGS_model)

# add pathway abbrev for heatmap
node_path_tbl = node_path_tbl %>%
  mutate(path_abbrev = case_when(
    stringr::str_detect(pathway, 'MAPK signaling') ~ 'MAPK',
    stringr::str_detect(pathway, 'Rho GTPases') ~ 'Rho',
    stringr::str_detect(pathway, 'Wnt signaling') ~ 'Wnt',
    stringr::str_detect(pathway, 'TGF-beta signaling') ~ 'TGF-b',
    stringr::str_detect(pathway, 'NF-kappa B signaling') ~ 'NF-kB',
    stringr::str_detect(pathway, 'mTOR signaling') ~ 'mTOR',
    stringr::str_detect(pathway, 'Apoptosis') ~ 'Apoptosis',
    stringr::str_detect(pathway, 'Cell cycle') ~ 'Cell Cycle',
    stringr::str_detect(pathway, 'Jak-STAT signaling') ~ 'JAK-STAT',
    stringr::str_detect(pathway, 'Receptor Tyrosine Kinases') ~ 'RTK',
    stringr::str_detect(pathway, 'PI3K-Akt signaling') ~ 'PI3K-AKT'))

saveRDS(object = node_path_tbl, file = 'data/node_path_tbl.rds')
