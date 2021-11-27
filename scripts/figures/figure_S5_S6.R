library(dplyr)
library(emba)
library(tibble)
library(ggplot2)
library(ggpubr)

###########################################
# CASCADE 1.0 topology scrambling results #
###########################################

# results from the scrambled topology simulations
# see 'scripts/get_syn_res_scrambled_topo_cascade1.R'
scrambled_topo_res = readRDS(file = 'data/scrambled_topo_res_cascade1.rds')

# results from bootstrap parameterization analysis
# see 'scripts/get_syn_res_boot_ss_cascade1.R'
boot_cascade1_res = readRDS(file = 'data/boot_cascade1_res.rds')

scrambled_topo_res = dplyr::bind_rows(scrambled_topo_res, boot_cascade1_res)

# group results by similarity score ('sim') and add group levels ('grp')
scrambled_topo_res =
  scrambled_topo_res %>% mutate(grp = factor(x =
      case_when(sim >= 0 & sim < 0.25 ~ '0 - 0.25',
        sim >= 0.25 & sim < 0.5 ~ '0.25 - 0.5',
        sim >= 0.5 & sim < 0.75 ~ '0.5 - 0.75',
        sim >= 0.75 & sim < 0.85 ~ '0.75 - 0.85',
        sim >= 0.85 & sim < 0.95 ~ '0.85 - 0.95',
        sim >= 0.95 & sim < 1 ~ '0.95 - 1',
        sim == 1 ~ 'Curated'),
    levels = c('0 - 0.25', '0.25 - 0.5', '0.5 - 0.75',
      '0.75 - 0.85', '0.85 - 0.95','0.95 - 1', 'Curated')))

##################################
# Figure S5 - ROC results in one #
##################################

# for reproducibility for `geom_jitter()`
set.seed(42)

## Source Scrambling
source_roc = scrambled_topo_res %>%
  filter(scramble_type == 'source' | scramble_type == 'none', !is.na(roc_auc)) %>%
  ggplot(aes(x = grp, y = roc_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,1)) +
  labs(x = 'Similarity Score to CASCADE 1.0 Topology', y = 'ROC AUC',
    title = "Source Scrambling vs Performance (ROC)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.7, y = 0.45, label = "Random (AUC = 0.5)")) +
  theme(plot.title = element_text(hjust = 0.5))

## Target Scrambling
target_roc = scrambled_topo_res %>%
  filter(scramble_type == 'target' | scramble_type == 'none', !is.na(roc_auc)) %>%
  ggplot(aes(x = grp, y = roc_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,1)) +
  labs(x = 'Similarity Score to CASCADE 1.0 Topology', y = 'ROC AUC',
    title = "Target Scrambling vs Performance (ROC)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.7, y = 0.45, label = "Random (AUC = 0.5)")) +
  theme(plot.title = element_text(hjust = 0.5))

## Sign Inversion
sign_roc = scrambled_topo_res %>%
  filter(scramble_type == 'effect' | scramble_type == 'none', !is.na(roc_auc)) %>%
  ggplot(aes(x = grp, y = roc_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,1)) +
  labs(x = 'Similarity Score to CASCADE 1.0 Topology', y = 'ROC AUC',
    title = "Sign Inversion vs Performance (ROC)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.7, y = 0.45, label = "Random (AUC = 0.5)")) +
  theme(plot.title = element_text(hjust = 0.5))

# no data points in the (0.95-1) class so we need to manually set the colors
# for consistency across all figures
set1_cols = RColorBrewer::brewer.pal(n = 7, name = 'Set1')[c(1:5,7)]

## All 3 types of topology scrambling
all_roc = scrambled_topo_res %>%
  filter(scramble_type == 'all' | scramble_type == 'none', !is.na(roc_auc)) %>%
  ggplot(aes(x = grp, y = roc_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = set1_cols) +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,1)) +
  labs(x = 'Similarity Score to CASCADE 1.0 Topology', y = 'ROC AUC',
    title = "All types of Scrambling vs Performance (ROC)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 5.8, y = 0.45, label = "Random (AUC = 0.5)")) +
  theme(plot.title = element_text(hjust = 0.5))

ggpubr::ggarrange(source_roc, target_roc, sign_roc, all_roc,
  labels = LETTERS[1:4], ncol = 2, nrow = 2)
ggplot2::ggsave(filename = 'scripts/figures/figure_S5_ROC.pdf', width = 14, height = 10, device = cairo_pdf)

#################################
# Figure S6 - PR results in one #
#################################

# Observed synergies (CASCADE 1.0)
ss_bliss_ew_file = 'results/link-only/cascade_1.0_ss_50sim_fixpoints_bliss_ensemblewise_synergies.tab'
ss_bliss_ew_synergies = emba::get_synergy_scores(ss_bliss_ew_file)
observed_synergies_file = 'data/observed_synergies_cascade_1.0'
observed_synergies = emba::get_observed_synergies(observed_synergies_file)
# 1 (positive/observed synergy) or 0 (negative/not observed) for all tested drug combinations
observed = sapply(ss_bliss_ew_synergies$perturbation %in% observed_synergies, as.integer)

# for reproducibility for `geom_jitter()`
set.seed(42)

## Source Scrambling
source_pr = scrambled_topo_res %>%
  filter(scramble_type == 'source' | scramble_type == 'none', !is.na(pr_auc)) %>%
  ggplot(aes(x = grp, y = pr_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,1)) +
  labs(x = 'Similarity Score to CASCADE 1.0 Topology', y = 'PR AUC',
    title = "Source Scrambling vs Performance (Precision-Recall)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = sum(observed)/length(observed), linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.7, y = 0.15, label = "Random (AUC = 0.2)")) +
  theme(plot.title = element_text(hjust = 0.5))

## Target Scrambling
target_pr = scrambled_topo_res %>%
  filter(scramble_type == 'target' | scramble_type == 'none', !is.na(pr_auc)) %>%
  ggplot(aes(x = grp, y = pr_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,1)) +
  labs(x = 'Similarity Score to CASCADE 1.0 Topology', y = 'PR AUC',
    title = "Target Scrambling vs Performance (Precision-Recall)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = sum(observed)/length(observed), linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.7, y = 0.13, label = "Random (AUC = 0.2)")) +
  theme(plot.title = element_text(hjust = 0.5))

## Sign Inversion
sign_pr = scrambled_topo_res %>%
  filter(scramble_type == 'effect' | scramble_type == 'none', !is.na(pr_auc)) %>%
  ggplot(aes(x = grp, y = pr_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,1)) +
  labs(x = 'Similarity Score to CASCADE 1.0 Topology', y = 'PR AUC',
    title = "Sign Inversion vs Performance (Precision-Recall)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = sum(observed)/length(observed), linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.9, y = 0.15, label = "Random (AUC = 0.2)"), size = 3) +
  theme(plot.title = element_text(hjust = 0.5))

## All 3 types of topology scrambling
all_pr = scrambled_topo_res %>%
  filter(scramble_type == 'all' | scramble_type == 'none', !is.na(pr_auc)) %>%
  ggplot(aes(x = grp, y = pr_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = set1_cols) +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,1)) +
  labs(x = 'Similarity Score to CASCADE 1.0 Topology', y = 'PR AUC',
    title = "All types of Scrambling vs Performance (Precision-Recall)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = sum(observed)/length(observed), linetype = 'dashed', color = "red") +
  geom_text(aes(x = 5.8, y = 0.15, label = "Random (AUC = 0.2)")) +
  theme(plot.title = element_text(hjust = 0.5))

ggpubr::ggarrange(source_pr, target_pr, sign_pr, all_pr,
  labels = LETTERS[1:4], ncol = 2, nrow = 2)
ggplot2::ggsave(filename = 'scripts/figures/figure_S6_PR.pdf', width = 14, height = 10, device = cairo_pdf)
