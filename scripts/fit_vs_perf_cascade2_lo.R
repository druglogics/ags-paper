###################################################
# Load Fitness vs Ensemble Model Performance Data #
# CASCADE 2.0, link-operator mutated models       #
###################################################
library(dplyr)
library(stringr)
library(emba)
library(tibble)
library(usefun)
library(PRROC)

# see Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `fit-vs-performance-results-bliss.tar.gz`
# get flipped training data results
data_dir = "/home/john/tmp/ags-paper/fit-vs-performance-results-bliss"

# get the AGS steady state (see `lo_mutated_models_heatmaps.R`)
steady_state = readRDS(file = 'data/steady_state.rds')

# define `beta` value for normalization
beta = -1

# get the prediction of the random proliferative models (150 simulations)
prolif_bliss_ensemblewise_150sim_file = paste0("results/link-only/cascade_2.0_rand_150sim_fixpoints_bliss_ensemblewise_synergies.tab")
prolif_bliss_ensemblewise_synergies_150sim = emba::get_synergy_scores(prolif_bliss_ensemblewise_150sim_file)

# Observed synergies for CASCADE 2.0
observed_synergies_file = 'data/observed_synergies_cascade_2.0'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)

data_list = list()
index = 1
for (res_dir in list.dirs(data_dir, recursive = FALSE)) {
  ew_synergies_file = list.files(path = res_dir, pattern = "ensemblewise_synergies", full.names = TRUE)
  ew_ss_scores = emba::get_synergy_scores(ew_synergies_file)

  # get number of nodes that were flipped in the steady state from directory name
  n_flips = stringr::str_match(res_dir, 'training_(.*)_flips')[1,2] %>% as.integer()

  # get the models stable states
  models_dir = paste0(res_dir, "/models")

  # you get messages for models with {#stable states} != 1
  # a few models have 2 stable states and they are not included
  # in the returned data frame
  models_stable_states = emba::get_stable_state_from_models_dir(models_dir)

  # calculate models fitness to AGS steady state
  models_fit = apply(models_stable_states[, names(steady_state)], 1,
    usefun::get_percentage_of_matches, steady_state)

  # calculate normalized model performance (ROC-AUC and PR-AUC)
  pred = dplyr::bind_cols(
    random = prolif_bliss_ensemblewise_synergies_150sim %>% select(score) %>% rename(random_score = score),
    ss = ew_ss_scores %>% select(score) %>% rename(ss_score = score),
    tibble::tibble(observed = sapply(ew_ss_scores$perturbation %in% observed_synergies, as.integer)))

  # get the normalized synergy scores
  pred = pred %>% mutate(combined_score = ss_score + beta * random_score)

  res_roc = PRROC::roc.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
    weights.class0 = pred %>% pull(observed))
  res_pr = PRROC::pr.curve(scores.class0 = pred %>% pull(combined_score) %>% (function(x) {-x}),
    weights.class0 = pred %>% pull(observed))

  # bind all to one (OneForAll)
  df = dplyr::bind_cols(
    roc_auc = res_roc$auc,
    pr_auc = res_pr$auc.davis.goadrich,
    avg_fit = mean(models_fit),
    per_flipped_data = n_flips/length(steady_state))

  print(paste('File', ew_synergies_file, 'fit-avg', mean(models_fit)))

  data_list[[index]] = df
  index = index + 1
}

res = dplyr::bind_rows(data_list)
saveRDS(res, file = "data/res_fit_aucs.rds")
