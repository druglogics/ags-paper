library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(emba)
library(PRROC)

# directory with the scrambled CASCADE 1.0 topology files
# see Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `scrambled_topologies_cascade_1.0.tar.gz`
topologies_dir = '/home/john/tmp/ags-paper/scrambled_topologies_cascade1'

# directory with the results of the script 'run_druglogics_synergy_scrambled_topo_cascade1.sh'
# See Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `synergy-res-scrambled-topo-cascade1.tar.gz`
syn_res_dir = '/home/john/tmp/ags-paper/synergy-res-scrambled-topo-cascade1'

# CASCADE 1.0 topology
edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_1.0.sif', delim = "\t", col_names = c('source', 'effect', 'target'), col_types = "ccc")

# observed synergies from (Flobak et al. 2015)
obs_syn_file = 'data/observed_synergies_cascade_1.0'
obs_syn = emba::get_observed_synergies(obs_syn_file)

# to measure similarity (reverse of scrambled-ness), compare edge full annotation
edge_annot = edge_tbl %>%
  mutate(edge = paste(source, effect, target)) %>% pull(edge)

data_list = list()
index = 1
for (topo_file in list.files(path = topologies_dir, full.names = TRUE)) {
  # get scrambled edges
  edge_tbl_scrambled = readr::read_delim(file = topo_file, delim = "\t",
    col_names = c('source', 'effect', 'target'), col_types = "ccc")
  edge_annot_scrambled = edge_tbl_scrambled %>%
    mutate(edge = paste(source, effect, target)) %>% pull(edge)

  # get similarity score
  sim = sum(edge_annot_scrambled %in% edge_annot)/length(edge_annot)

  # get topology pattern to search in the synergy results
  topo_pat = stringr::str_extract(string = topo_file, pattern = "(cascade1_.*.sif)")
  topo_pat = stringr::str_remove(string = topo_pat, pattern = ".sif")
  topo_pat = stringr::str_c(topo_pat, "_", sep = "")

  # get type of scrambling (source node, target node, effect reversal, or all together?)
  if (stringr::str_detect(string = topo_pat, pattern = "src_permut"))
    type = "source"
  else if (stringr::str_detect(string = topo_pat, pattern = "trg_permut"))
    type = "target"
  else if (stringr::str_detect(string = topo_pat, pattern = "eff_permut"))
    type = "effect"
  else if (stringr::str_detect(string = topo_pat, pattern = "all_permut"))
    type = "all"

  res_dirs = list.files(path = syn_res_dir, pattern = topo_pat, full.names = TRUE)

  # check: one dir with calibrated to steady state results, one with random proliferative results
  stopifnot(length(res_dirs) == 2)
  random_dir = res_dirs[stringr::str_which(string = res_dirs, pattern = "random")]
  ss_dir     = res_dirs[stringr::str_which(string = res_dirs, pattern = "ss")]

  # ensemble-wise synergy files
  ew_random_file = list.files(path = random_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
  ew_ss_file     = list.files(path = ss_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
  if (length(ew_random_file) == 0 | length(ew_ss_file) == 0) {
    df = dplyr::bind_cols(sim = sim, scramble_type = type, roc_auc = NA, pr_auc = NA)
    data_list[[index]] = df
    index = index + 1

    next
  }

  # get ensemble-wise synergy scores
  random_scores = emba::get_synergy_scores(file_name = ew_random_file)
  ss_scores     = emba::get_synergy_scores(file_name = ew_ss_file)

  # get observed synergy vector (1 = synergy, 0 = no synergy)
  stopifnot(all(ss_scores$perturbation == random_scores$perturbation))
  observed = sapply(ss_scores$perturbation %in% obs_syn, as.integer)

  # tidy data
  pred = dplyr::bind_cols(
    random = random_scores %>% select(perturbation, score) %>% rename(random_score = score),
    ss = ss_scores %>% select(score) %>% rename(ss_score = score),
    tibble::as_tibble_col(observed, column_name = "observed"))

  # calculate normalized synergy scores
  pred = pred %>% mutate(combined_score = ss_score - random_score)

  # get ROC and PR AUC results
  res_roc = PRROC::roc.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
    weights.class0 = pred %>% pull(observed))
  res_pr = PRROC::pr.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
    weights.class0 = pred %>% pull(observed))

  # bind all to one (OneForAll)
  df = dplyr::bind_cols(sim = sim, scramble_type = type,
    roc_auc = res_roc$auc, pr_auc = res_pr$auc.davis.goadrich)

  data_list[[index]] = df
  index = index + 1
}

scrambled_topo_res_cascade1 = dplyr::bind_rows(data_list)
saveRDS(object = scrambled_topo_res_cascade1, file = 'data/scrambled_topo_res_cascade1.rds')
