###########################################################
# Compare performance from two link-operator model pools: #
# active vs inhibited ERK_f                               #
###########################################################
library(dplyr)
library(stringr)
library(emba)
library(tibble)
library(usefun)
library(PRROC)

# from http://tiny.cc/ags-paper-zenodo, file `erk_perf_investigation.tar.gz`,
data_dir = "/home/john/tmp/ags-paper/erk_perf_investigation"

# define `beta` value for normalization
beta = -1

# Observed synergies for CASCADE 2.0
observed_synergies_file = 'data/observed_synergies_cascade_2.0'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)

# random model predictions
# from http://tiny.cc/ags-paper-zenodo, file `parameterization-comp.tar.gz` directory `link-only`
random_ew_scores_file = '/home/john/tmp/ags-paper/parameterization-comp/link-only/cascade_2.0_rand_150sim_fixpoints_bliss_20200505_063817/cascade_2.0_rand_150sim_fixpoints_bliss_ensemblewise_synergies.tab'
random_ew_scores = emba::get_synergy_scores(random_ew_scores_file)

# define data list that is going to store all results
data_list = list()
index = 1

# drabme results from the 2 pools
for (res_dir in list.files(data_dir, recursive = FALSE, full.names = TRUE, pattern = 'cascade')) {
  # get ensemble-wise synergies scores
  ew_synergies_file = list.files(path = res_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
  ss_ew_scores = emba::get_synergy_scores(ew_synergies_file)

  # calculate normalized model performance (ROC-AUC and PR-AUC)
  pred = dplyr::bind_cols(
    random_ew_scores %>% rename(random_score = score),
    ss_ew_scores %>% select(score) %>% rename(ss_score = score),
    tibble::tibble(observed = sapply(ss_ew_scores$perturbation %in% observed_synergies, as.integer)))

  # get the normalized synergy scores
  pred = pred %>% mutate(combined_score = ss_score + beta * random_score)

  print(ifelse(stringr::str_detect(string = res_dir, pattern = 'active'), 'active', 'inhibited'))
  print(pred %>% filter(observed == 1) %>% arrange(combined_score))

  res_roc = PRROC::roc.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
    weights.class0 = pred %>% pull(observed))
  res_pr = PRROC::pr.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
    weights.class0 = pred %>% pull(observed))

  # bind all to one (OneForAll)
  df = dplyr::bind_cols(roc_auc = res_roc$auc, pr_auc = res_pr$auc.davis.goadrich,
    erk_state = ifelse(stringr::str_detect(string = res_dir, pattern = 'active'), 'active', 'inhibited'))

  data_list[[index]] = df
  index = index + 1
}

res = dplyr::bind_rows(data_list)
saveRDS(res, file = "data/res_erk.rds")
