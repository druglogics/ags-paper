library(emba)
library(usefun)
library(dplyr)
library(tibble)
library(PRROC)

# Load synegy predictions (CASCADE 1.0, 50 simulations, Bliss)
# 'ss' => calibrated models, 'prolif' => random models
ss_bliss_ew_file     = 'results/link-only/cascade_1.0_ss_50sim_fixpoints_bliss_ensemblewise_synergies.tab'
prolif_bliss_ew_file = 'results/link-only/cascade_1.0_rand_50sim_fixpoints_bliss_ensemblewise_synergies.tab'

ss_bliss_ew_synergies     = emba::get_synergy_scores(ss_bliss_ew_file)
prolif_bliss_ew_synergies = emba::get_synergy_scores(prolif_bliss_ew_file)

all(ss_bliss_ew_synergies$perturbation == prolif_bliss_ew_synergies$perturbation)

# Observed synergies (CASCADE 1.0)
observed_synergies_file = 'data/observed_synergies_cascade_1.0'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)
# 1 (positive/observed synergy) or 0 (negative/not observed) for all tested drug combinations
observed = sapply(ss_bliss_ew_synergies$perturbation %in% observed_synergies, as.integer)

# make predictions table and save to file
beta = -1
pred_ew_bliss = dplyr::bind_cols(
  ss_bliss_ew_synergies %>% rename(ss_score = score),
  prolif_bliss_ew_synergies %>% select(score) %>% rename(prolif_score = score),
  tibble::as_tibble_col(observed, column_name = 'observed')) %>%
  mutate(combined_score = ss_score + beta * prolif_score, .before = observed)
readr::write_csv(x = pred_ew_bliss, file = 'scripts/figures/cascade1.0_prediction_results.csv')

# ROC statistics
res_ss_ew = usefun::get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score", label_col = "observed")
res_prolif_ew = usefun::get_roc_stats(df = pred_ew_bliss, pred_col = "prolif_score", label_col = "observed")
res_comb_pred = usefun::get_roc_stats(df = pred_ew_bliss, pred_col = "combined_score", label_col = "observed")

# Table to view the combined ROC statistics
DT::datatable(data = res_comb_pred$roc_stats, options =
    list(pageLength = 5, lengthMenu = c(5, 10, 16), searching = FALSE)) %>%
  DT::formatRound(c(1,6,7,8,9), digits = 3)

set1_cols = RColorBrewer::brewer.pal(n = 9, name = "Set1")

# Figure 2C - ROC (1)
# 'Calibrated' here is Calibrated + normalized to random proliferative model predictions
# i.e. the combined predictor
png(filename = 'scripts/figures/figure_2C_1.png', width = 5, height = 5, units = 'in', res = 300)
plot(x = res_comb_pred$roc_stats$FPR, y = res_comb_pred$roc_stats$TPR,
  type = 'l', lwd = 3, col = set1_cols[1], main = 'ROC curve, Ensemble-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_prolif_ew$roc_stats$FPR, y = res_prolif_ew$roc_stats$TPR,
  lwd = 3, col = set1_cols[2])
legend('bottomright', title = 'AUC', col = set1_cols[1:2], pch = 19,
  legend = c(paste(round(res_comb_pred$AUC, digits = 2), "Calibrated"),
    paste(round(res_prolif_ew$AUC, digits = 2), "Random")), cex = 1.3)
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
dev.off()

# PR statistics
res_ss_ew_pr     = PRROC::pr.curve(scores.class0 = pred_ew_bliss %>%
    pull(ss_score) %>% (function(x) {-x}),
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)
res_prolif_ew_pr = PRROC::pr.curve(scores.class0 = pred_ew_bliss %>%
    pull(prolif_score) %>% (function(x) {-x}),
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)
res_comb_pred_pr = PRROC::pr.curve(scores.class0 = pred_ew_bliss %>%
    pull(combined_score) %>% (function(x) {-x}),
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)

# Figure 2C - PR (2)
# 'Calibrated' here is Calibrated + normalized to random proliferative model predictions
# i.e. the combined predictor
png(filename = 'scripts/figures/figure_2C_2.png', width = 5, height = 5, units = 'in', res = 300)
plot(res_comb_pred_pr, main = 'PR curve, Ensemble-wise synergies (Bliss)',
  auc.main = FALSE, color = set1_cols[1], rand.plot = TRUE)
plot(res_prolif_ew_pr, add = TRUE, color = set1_cols[2])
legend(x = 0, y = 0.9, title = 'AUC', col = set1_cols[1:2], pch = 19, cex = 1.3,
  legend = c(paste(round(res_comb_pred_pr$auc.davis.goadrich, digits = 2), "Calibrated"),
    paste(round(res_prolif_ew_pr$auc.davis.goadrich, digits = 2), "Random")))
grid(lwd = 0.5)
dev.off()
