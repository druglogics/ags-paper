library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(emba)
library(PRROC)

# see Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `ss_cascade1_model_bootstrap.tar.gz`
data_dir = '/home/john/tmp/ags-paper/ss_cascade1_model_bootstrap'
res_dirs = list.files(path = data_dir, pattern = 'bliss', full.names = TRUE)
stopifnot(length(res_dirs) == 52)

# get random proliferative results (50 simulations)
random_dir = res_dirs[stringr::str_which(string = res_dirs, pattern = "rand")]
ew_random_file = list.files(path = random_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
random_scores = emba::get_synergy_scores(file_name = ew_random_file)

# observed synergies from (Flobak et al. 2015)
obs_syn_file = 'data/observed_synergies_cascade_1.0'
obs_syn = emba::get_observed_synergies(obs_syn_file)
observed = sapply(random_scores$perturbation %in% obs_syn, as.integer)

# keep only the results from the bootstrap drabme analysis (`boostrap_models_drabme_cascade1.sh`)
ss_dirs = res_dirs[stringr::str_detect(string = res_dirs, pattern = 'bliss_batch')]
stopifnot(length(ss_dirs) == 50) # 50 batches

data_list = list()
index = 1
for (ss_dir in ss_dirs) {
  # get the calibrated to steady state result from the specific batch
  ew_ss_file = list.files(path = ss_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
  ss_scores = emba::get_synergy_scores(file_name = ew_ss_file)

  # data check
  stopifnot(random_scores$perturbation == ss_scores$perturbation)

  # tidy data
  pred = dplyr::bind_cols(
    random_scores %>% select(perturbation, score) %>% rename(random_score = score),
    ss_scores %>% select(score) %>% rename(ss_score = score))

  # calculate normalized synergy scores
  pred = pred %>% mutate(combined_score = ss_score - random_score)

  # get ROC and PR AUC results
  res_roc = PRROC::roc.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
    weights.class0 = observed)
  res_pr = PRROC::pr.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
    weights.class0 = observed)

  # bind all to one (OneForAll)
  # `sim` = 1 because all these results come from the original,
  # curated CASCADE 1.0 topology, so `scramble_type` is 'none'
  df = dplyr::bind_cols(sim = 1, scramble_type = 'none',
    roc_auc = res_roc$auc, pr_auc = res_pr$auc.davis.goadrich)

  data_list[[index]] = df
  index = index + 1
}

boot_cascade1_res = dplyr::bind_rows(data_list)

saveRDS(object = boot_cascade1_res, file = 'data/boot_cascade1_res.rds')
