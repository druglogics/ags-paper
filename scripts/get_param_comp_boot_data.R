##############################################################
# Load Bootstrap Models Parameterization vs Performance Data #
# CASCADE 2.0, 3 parameterization schemes used               #
##############################################################
library(dplyr)
library(stringr)
library(emba)
library(tibble)
library(usefun)
library(PRROC)

# from http://tiny.cc/ags-paper-zenodo, file `parameterization-comp.tar.gz`
# lo = 'link-operator mutations', topo = 'topology mutations', both = 'link-operator and topology mutations'
data_dir_lo = "/home/john/tmp/ags-paper/parameterization-comp/link-only/"
data_dir_topo = "/home/john/tmp/ags-paper/parameterization-comp/topology-only/"
data_dir_both = "/home/john/tmp/ags-paper/parameterization-comp/topo-and-link/"

# define `beta` value for normalization
beta = -1

# Observed synergies for CASCADE 2.0
observed_synergies_file = 'data/observed_synergies_cascade_2.0'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)

# define data list that is going to store all results
data_list = list()
index = 1

## Link-operator only Mutations

# random model predictions
random_ew_scores_file_lo = paste0(data_dir_lo, "cascade_2.0_rand_150sim_fixpoints_bliss_20200505_063817/cascade_2.0_rand_150sim_fixpoints_bliss_ensemblewise_synergies.tab")
random_ew_scores_lo = emba::get_synergy_scores(random_ew_scores_file_lo)

# calibrated bootstrap data
for (res_dir in list.dirs(paste0(data_dir_lo, "boot_res"), recursive = FALSE)) {
  # check only the simulation (not the `models_batch_*`) directories
  if (stringr::str_detect(string = res_dir, pattern = "cascade_2.0_ss_bliss_batch")) {
    ew_synergies_file = list.files(path = res_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
    ss_ew_scores = emba::get_synergy_scores(ew_synergies_file)

    # calculate normalized model performance (ROC-AUC and PR-AUC)
    pred = dplyr::bind_cols(
      random_ew_scores_lo %>% select(score) %>% rename(random_score = score),
      ss_ew_scores %>% select(score) %>% rename(ss_score = score),
      tibble::tibble(observed = sapply(ss_ew_scores$perturbation %in% observed_synergies, as.integer)))

    # get the normalized synergy scores
    pred = pred %>% mutate(combined_score = ss_score + beta * random_score)

    res_roc = PRROC::roc.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
      weights.class0 = pred %>% pull(observed))
    res_pr = PRROC::pr.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
      weights.class0 = pred %>% pull(observed))

    # bind all to one (OneForAll)
    df = dplyr::bind_cols(roc_auc = res_roc$auc, pr_auc = res_pr$auc.davis.goadrich,
      param = "link-only")

    data_list[[index]] = df
    index = index + 1
  }
}

## Topology-only mutations

# random model predictions
random_ew_scores_file_topo = paste0(data_dir_topo, "cascade_2.0_rand_150sim_fixpoints_bliss_20200429_023822/cascade_2.0_rand_150sim_fixpoints_bliss_ensemblewise_synergies.tab")
random_ew_scores_topo = emba::get_synergy_scores(random_ew_scores_file_topo)

# calibrated bootstrap data
for (res_dir in list.dirs(paste0(data_dir_topo, "boot_res"), recursive = FALSE)) {
  # check only the simulation (not the `models_batch_*`) directories
  if (stringr::str_detect(string = res_dir, pattern = "cascade_2.0_ss_bliss_batch")) {
    ew_synergies_file = list.files(path = res_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
    ss_ew_scores = emba::get_synergy_scores(ew_synergies_file)

    # calculate normalized model performance (ROC-AUC and PR-AUC)
    pred = dplyr::bind_cols(
      random_ew_scores_topo %>% select(score) %>% rename(random_score = score),
      ss_ew_scores %>% select(score) %>% rename(ss_score = score),
      tibble::tibble(observed = sapply(ss_ew_scores$perturbation %in% observed_synergies, as.integer)))

    # get the normalized synergy scores
    pred = pred %>% mutate(combined_score = ss_score + beta * random_score)

    res_roc = PRROC::roc.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
      weights.class0 = pred %>% pull(observed))
    res_pr = PRROC::pr.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
      weights.class0 = pred %>% pull(observed))

    # bind all to one (OneForAll)
    df = dplyr::bind_cols(roc_auc = res_roc$auc, pr_auc = res_pr$auc.davis.goadrich,
      param = "topology-only")

    data_list[[index]] = df
    index = index + 1
  }
}

## Both Link-operator and Topology mutations

# random model predictions
random_ew_scores_file_both = paste0(data_dir_both, "cascade_2.0_rand_150sim_fixpoints_bliss_20200430_122450/cascade_2.0_rand_150sim_fixpoints_bliss_ensemblewise_synergies.tab")
random_ew_scores_both = emba::get_synergy_scores(random_ew_scores_file_both)

# calibrated bootstrap data
for (res_dir in list.dirs(paste0(data_dir_both, "boot_res"), recursive = FALSE)) {
  # check only the simulation (not the `models_batch_*`) directories
  if (stringr::str_detect(string = res_dir, pattern = "cascade_2.0_ss_bliss_batch")) {
    ew_synergies_file = list.files(path = res_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
    ss_ew_scores = emba::get_synergy_scores(ew_synergies_file)

    # calculate normalized model performance (ROC-AUC and PR-AUC)
    pred = dplyr::bind_cols(
      random_ew_scores_both %>% select(score) %>% rename(random_score = score),
      ss_ew_scores %>% select(score) %>% rename(ss_score = score),
      tibble::tibble(observed = sapply(ss_ew_scores$perturbation %in% observed_synergies, as.integer)))

    # get the normalized synergy scores
    pred = pred %>% mutate(combined_score = ss_score + beta * random_score)

    res_roc = PRROC::roc.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
      weights.class0 = pred %>% pull(observed))
    res_pr = PRROC::pr.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
      weights.class0 = pred %>% pull(observed))

    # bind all to one (OneForAll)
    df = dplyr::bind_cols(roc_auc = res_roc$auc, pr_auc = res_pr$auc.davis.goadrich,
      param = "topo-and-link")

    data_list[[index]] = df
    index = index + 1
  }
}

res = dplyr::bind_rows(data_list)
saveRDS(res, file = "data/res_param_boot_aucs.rds")
