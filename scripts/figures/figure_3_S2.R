library(emba)
library(usefun)
library(dplyr)
library(tibble)
library(PRROC)
library(pROC)
source('scripts/figures/helpers.R')

# Load synegy predictions (CASCADE 2.0, 150 simulations, Bliss)
# 'ss' => calibrated models, 'prolif' => random models
ss_bliss_ew_file     = 'results/link-only/cascade_2.0_ss_150sim_fixpoints_bliss_ensemblewise_synergies.tab'
prolif_bliss_ew_file = 'results/link-only/cascade_2.0_rand_150sim_fixpoints_bliss_ensemblewise_synergies.tab'

ss_bliss_ew_synergies     = emba::get_synergy_scores(ss_bliss_ew_file)
prolif_bliss_ew_synergies = emba::get_synergy_scores(prolif_bliss_ew_file)

all(ss_bliss_ew_synergies$perturbation == prolif_bliss_ew_synergies$perturbation)

# Observed synergies (CASCADE 2.0)
observed_synergies_file = 'data/observed_synergies_cascade_2.0'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)
# 1 (positive/observed synergy) or 0 (negative/not observed) for all tested drug combinations
observed = sapply(ss_bliss_ew_synergies$perturbation %in% observed_synergies, as.integer)

# make predictions table and save to file
beta = -1
pred_ew_bliss = dplyr::bind_cols(
  ss_bliss_ew_synergies %>% rename(ss_score = score), # calibrated models synergy scores
  prolif_bliss_ew_synergies %>% select(score) %>% rename(prolif_score = score),  # random proliferative models synergy scores
  tibble::as_tibble_col(observed, column_name = 'observed')) %>%
  mutate(combined_score = ss_score + beta * prolif_score, .before = observed) # calibrated normalized synergy scores
readr::write_csv(x = pred_ew_bliss, file = 'scripts/figures/cascade2.0_prediction_results.csv')

# ROC statistics
res_ss_ew = usefun::get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score", label_col = "observed")
res_prolif_ew = usefun::get_roc_stats(df = pred_ew_bliss, pred_col = "prolif_score", label_col = "observed")
res_comb_pred = usefun::get_roc_stats(df = pred_ew_bliss, pred_col = "combined_score", label_col = "observed")

# Table to view the combined ROC statistics
DT::datatable(data = res_comb_pred$roc_stats, options =
    list(pageLength = 5, lengthMenu = c(5, 10, 16), searching = FALSE)) %>%
  DT::formatRound(c(1,6,7,8,9), digits = 5)

# Compare ROC AUCs
roc_prolif = pROC::roc(response = pred_ew_bliss$observed,
  predictor = pred_ew_bliss$prolif_score, direction = ">")
# controls > cases (lower more synergistic)
roc_combined = pROC::roc(response = pred_ew_bliss$observed,
  predictor = pred_ew_bliss$combined_score, direction = ">")
set.seed(42)
pROC::roc.test(roc_prolif, roc_combined, method = "bootstrap", boot.n = 10000,
  direction = ">")

set1_cols = RColorBrewer::brewer.pal(n = 9, name = "Set1")

# Figure 3 - ROC
# 'Calibrated' here is Calibrated + normalized to random proliferative model predictions
# i.e. the combined predictor
cairo_pdf(filename = 'scripts/figures/figure_3_ROC.pdf', width = 5, height = 5)
plot(x = res_comb_pred$roc_stats$FPR, y = res_comb_pred$roc_stats$TPR,
  type = 'l', lwd = 3, col = set1_cols[1], main = 'ROC curve, Ensemble-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_prolif_ew$roc_stats$FPR, y = res_prolif_ew$roc_stats$TPR,
  lwd = 3, col = set1_cols[2])
legend('topleft', title = 'AUC', col = set1_cols[1:2], pch = 19,
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

# Compare PR AUCs
set.seed(42)
pr.test(labels = pred_ew_bliss$observed, pred1 = -pred_ew_bliss$prolif_score,
  pred2 = -pred_ew_bliss$combined_score, boot.n = 10000)

# Figure 3 - PR
# 'Calibrated' here is Calibrated + normalized to random proliferative model predictions
# i.e. the combined predictor
cairo_pdf(filename = 'scripts/figures/figure_3_PR.pdf', width = 5, height = 5)
plot(res_comb_pred_pr, main = 'PR curve, Ensemble-wise synergies (Bliss)',
  auc.main = FALSE, color = set1_cols[1], rand.plot = TRUE)
plot(res_prolif_ew_pr, add = TRUE, color = set1_cols[2])
legend(x = 0, y = 0.9, title = 'AUC', col = set1_cols[1:2], pch = 19, cex = 1.3,
  legend = c(paste(round(res_comb_pred_pr$auc.davis.goadrich, digits = 2), "Calibrated"),
    paste(round(res_prolif_ew_pr$auc.davis.goadrich, digits = 2), "Random")))
grid(lwd = 0.5)
dev.off()

# Figure 3 - ROC and PR combined
cairo_pdf(filename = 'scripts/figures/figure_3.pdf', width = 10, height = 5)
par(mfrow = c(1,2))
plot(x = res_comb_pred$roc_stats$FPR, y = res_comb_pred$roc_stats$TPR,
  type = 'l', lwd = 3, col = set1_cols[1], main = 'ROC curve, Ensemble-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_prolif_ew$roc_stats$FPR, y = res_prolif_ew$roc_stats$TPR,
  lwd = 3, col = set1_cols[2])
legend('topleft', title = 'AUC', col = set1_cols[1:2], pch = 19,
  legend = c(paste(round(res_comb_pred$AUC, digits = 2), "Calibrated"),
    paste(round(res_prolif_ew$AUC, digits = 2), "Random")), cex = 1.3)
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(res_comb_pred_pr, main = 'PR curve, Ensemble-wise synergies (Bliss)',
  auc.main = FALSE, color = set1_cols[1], rand.plot = TRUE)
plot(res_prolif_ew_pr, add = TRUE, color = set1_cols[2])
legend(x = 0, y = 0.9, title = 'AUC', col = set1_cols[1:2], pch = 19, cex = 1.3,
  legend = c(paste(round(res_comb_pred_pr$auc.davis.goadrich, digits = 2), "Calibrated"),
    paste(round(res_prolif_ew_pr$auc.davis.goadrich, digits = 2), "Random")))
grid(lwd = 0.5)
dev.off()

# Figure S2 - ROC
# Calibrated non-normalized vs calibrated normalized vs random proliferative
cairo_pdf(filename = 'scripts/figures/figure_S2_ROC.pdf', width = 5, height = 5)
plot(x = res_comb_pred$roc_stats$FPR, y = res_comb_pred$roc_stats$TPR,
  type = 'l', lwd = 3, col = set1_cols[1],
  main = ('ROC curve, Ensemble-wise synergies (Bliss)'),
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_prolif_ew$roc_stats$FPR, y = res_prolif_ew$roc_stats$TPR,
  lwd = 3, col = set1_cols[2])
lines(x = res_ss_ew$roc_stats$FPR, y = res_ss_ew$roc_stats$TPR,
  lwd = 3, col = set1_cols[3])
legend('topleft', title = 'AUC', col = set1_cols[c(3,1,2)], pch = 19, cex = 1,
  legend = c(paste(round(res_ss_ew$AUC, digits = 2), 'Calibrated (non-normalized)'),
    paste(round(res_comb_pred$AUC, digits = 2), 'Calibrated'),
    paste(round(res_prolif_ew$AUC, digits = 2), 'Random')))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
dev.off()

# Figure S2 - PR
# Calibrated non-normalized vs calibrated normalized vs random proliferative
cairo_pdf(filename = 'scripts/figures/figure_S2_PR.pdf', width = 5, height = 5)
plot(res_comb_pred_pr, main = 'PR curve, Ensemble-wise synergies (Bliss)',
  auc.main = FALSE, color = set1_cols[1], rand.plot = TRUE, lwd = 3)
plot(res_prolif_ew_pr, add = TRUE, color = set1_cols[2], lwd = 3)
plot(res_ss_ew_pr, add = TRUE, color = set1_cols[3], lwd = 2)
legend('topright', title = 'AUC', col = set1_cols[c(3,1,2)], pch = 19, cex = 1,
  legend = c(paste(round(res_ss_ew_pr$auc.davis.goadrich, digits = 2), 'Calibrated (non-normalized)'),
    paste(round(res_comb_pred_pr$auc.davis.goadrich, digits = 2), 'Calibrated'),
    paste(round(res_prolif_ew_pr$auc.davis.goadrich, digits = 2), 'Random')))
grid(lwd = 0.5)
dev.off()

# Figure S2 - ROC and PR combined
cairo_pdf(filename = 'scripts/figures/figure_S2.pdf', width = 10, height = 5)
par(mfrow = c(1,2))
plot(x = res_comb_pred$roc_stats$FPR, y = res_comb_pred$roc_stats$TPR,
  type = 'l', lwd = 3, col = set1_cols[1],
  main = ('ROC curve, Ensemble-wise synergies (Bliss)'),
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_prolif_ew$roc_stats$FPR, y = res_prolif_ew$roc_stats$TPR,
  lwd = 3, col = set1_cols[2])
lines(x = res_ss_ew$roc_stats$FPR, y = res_ss_ew$roc_stats$TPR,
  lwd = 3, col = set1_cols[3])
legend('topleft', title = 'AUC', col = set1_cols[c(3,1,2)], pch = 19, cex = 1,
  legend = c(paste(round(res_ss_ew$AUC, digits = 2), 'Calibrated (non-normalized)'),
    paste(round(res_comb_pred$AUC, digits = 2), 'Calibrated'),
    paste(round(res_prolif_ew$AUC, digits = 2), 'Random')))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(res_comb_pred_pr, main = 'PR curve, Ensemble-wise synergies (Bliss)',
  auc.main = FALSE, color = set1_cols[1], rand.plot = TRUE, lwd = 3)
plot(res_prolif_ew_pr, add = TRUE, color = set1_cols[2], lwd = 3)
plot(res_ss_ew_pr, add = TRUE, color = set1_cols[3], lwd = 2)
legend('topright', title = 'AUC', col = set1_cols[c(3,1,2)], pch = 19, cex = 1,
  legend = c(paste(round(res_ss_ew_pr$auc.davis.goadrich, digits = 2), 'Calibrated (non-normalized)'),
    paste(round(res_comb_pred_pr$auc.davis.goadrich, digits = 2), 'Calibrated'),
    paste(round(res_prolif_ew_pr$auc.davis.goadrich, digits = 2), 'Random')))
grid(lwd = 0.5)
dev.off()
