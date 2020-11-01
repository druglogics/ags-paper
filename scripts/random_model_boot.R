#######################################
# Load bootstrap random model results #
#######################################
library(dplyr)
library(stringr)
library(emba)
library(PRROC)

# see Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `random_model_bootstrap.tar.gz`
data_dir = "/home/john/tmp/ags-paper/random_model_bootstrap"

# Observed synergies for CASCADE 2.0
observed_synergies_file = 'data/observed_synergies_cascade_2.0'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)

data_list = list()
index = 11
for (res_dir in list.dirs(data_dir, recursive = FALSE)) {
  if (stringr::str_detect(string = res_dir, pattern = "cascade_2.0_rand_prolif_bliss_batch")) {
    ew_synergies_file = list.files(path = res_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
    rand_scores = emba::get_synergy_scores(ew_synergies_file)
    observed = sapply(rand_scores$perturbation %in% observed_synergies, as.integer)

    res_roc = PRROC::roc.curve(scores.class0 = rand_scores$score %>% (function(x) {-x}),
      weights.class0 = observed)
    res_pr = PRROC::pr.curve(scores.class0 = rand_scores$score %>% (function(x) {-x}),
      weights.class0 = observed)

    # bind all to one (OneForAll)
    df = dplyr::bind_cols(roc_auc = res_roc$auc, pr_auc = res_pr$auc.davis.goadrich)
    data_list[[index]] = df
    index = index + 1
  }
}

rand_res = dplyr::bind_rows(data_list)
saveRDS(rand_res, file = "data/bootstrap_rand_res.rds")
