library(dplyr)
library(emba)
library(tibble)
library(ggplot2)
library(ggpubr)

###########################################
# CASCADE 2.0 topology scrambling results #
###########################################

# results from the scrambled topology simulations
# see 'scripts/get_syn_res_scrambled_topo_cascade2.R'
scrambled_topo_res = readRDS(file = 'data/scrambled_topo_res_cascade2.rds')

# results from bootstrap parameterization analysis
# see 'scripts/get_param_comp_boot_data.R'
boot_res = readRDS(file = "data/res_param_boot_aucs.rds")

# the un-scrambled topology results have a similarity score ('sim') equal to 1,
# and 'none' scrambling whatsoever as `scramble_type`
lo_boot_res = boot_res %>%
  tibble::add_column(sim = 1, scramble_type = 'none', .before = 1) %>%
  filter(param == 'link-only') %>% # keep only the link-operator results
  select(-one_of("param")) # remove unwanted column

scrambled_topo_res = dplyr::bind_rows(scrambled_topo_res, lo_boot_res)

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

#################################
# Figure 5 - ROC results in one #
#################################

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
  labs(x = 'Similarity Score to CASCADE 2.0 Topology', y = 'ROC AUC',
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
  labs(x = 'Similarity Score to CASCADE 2.0 Topology', y = 'ROC AUC',
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
  labs(x = 'Similarity Score to CASCADE 2.0 Topology', y = 'ROC AUC',
    title = "Sign Inversion vs Performance (ROC)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.7, y = 0.45, label = "Random (AUC = 0.5)")) +
  theme(plot.title = element_text(hjust = 0.5))

## All 3 types of topology scrambling
all_roc = scrambled_topo_res %>%
  filter(scramble_type == 'all' | scramble_type == 'none', !is.na(roc_auc)) %>%
  ggplot(aes(x = grp, y = roc_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,1)) +
  labs(x = 'Similarity Score to CASCADE 2.0 Topology', y = 'ROC AUC',
    title = "All types of Scrambling vs Performance (ROC)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.5, y = 0.45, label = "Random (AUC = 0.5)")) +
  theme(plot.title = element_text(hjust = 0.5))

ggpubr::ggarrange(source_roc, target_roc, sign_roc, all_roc,
  labels = LETTERS[1:4], ncol = 2, nrow = 2)
ggplot2::ggsave(filename = 'scripts/figures/figure_5_ROC.pdf', width = 14, height = 10, device = cairo_pdf)

#################################
# Figure S4 - PR results in one #
#################################

# Observed synergies (CASCADE 2.0)
ss_bliss_ew_file = 'results/link-only/cascade_2.0_ss_150sim_fixpoints_bliss_ensemblewise_synergies.tab'
ss_bliss_ew_synergies = emba::get_synergy_scores(ss_bliss_ew_file)
observed_synergies_file = 'data/observed_synergies_cascade_2.0'
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
  ylim(c(0,0.6)) +
  labs(x = 'Similarity Score to CASCADE 2.0 Topology', y = 'PR AUC',
    title = "Source Scrambling vs Performance (Precision-Recall)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = sum(observed)/length(observed), linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.7, y = 0.01, label = "Random (AUC = 0.04)")) +
  theme(plot.title = element_text(hjust = 0.5))

## Target Scrambling
target_pr = scrambled_topo_res %>%
  filter(scramble_type == 'target' | scramble_type == 'none', !is.na(pr_auc)) %>%
  ggplot(aes(x = grp, y = pr_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,0.6)) +
  labs(x = 'Similarity Score to CASCADE 2.0 Topology', y = 'PR AUC',
    title = "Target Scrambling vs Performance (Precision-Recall)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = sum(observed)/length(observed), linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.7, y = 0.01, label = "Random (AUC = 0.04)")) +
  theme(plot.title = element_text(hjust = 0.5))

## Sign Inversion
sign_pr = scrambled_topo_res %>%
  filter(scramble_type == 'effect' | scramble_type == 'none', !is.na(pr_auc)) %>%
  ggplot(aes(x = grp, y = pr_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,0.6)) +
  labs(x = 'Similarity Score to CASCADE 2.0 Topology', y = 'PR AUC',
    title = "Sign Inversion vs Performance (Precision-Recall)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = sum(observed)/length(observed), linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.5, y = 0.01, label = "Random (AUC = 0.04)")) +
  theme(plot.title = element_text(hjust = 0.5))

## All 3 types of topology scrambling
all_pr = scrambled_topo_res %>%
  filter(scramble_type == 'all' | scramble_type == 'none', !is.na(pr_auc)) %>%
  ggplot(aes(x = grp, y = pr_auc, fill = grp)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_brewer(palette = 'Set1') +
  geom_jitter(shape = 20, position = position_jitter(0.2), show.legend = FALSE) +
  ylim(c(0,0.6)) +
  labs(x = 'Similarity Score to CASCADE 2.0 Topology', y = 'PR AUC',
    title = "All types of Scrambling vs Performance (Precision-Recall)") +
  theme_classic(base_size = 14) +
  geom_hline(yintercept = sum(observed)/length(observed), linetype = 'dashed', color = "red") +
  geom_text(aes(x = 6.5, y = 0.01, label = "Random (AUC = 0.04)")) +
  theme(plot.title = element_text(hjust = 0.5))

ggpubr::ggarrange(source_pr, target_pr, sign_pr, all_pr,
  labels = LETTERS[1:4], ncol = 2, nrow = 2)
ggplot2::ggsave(filename = 'scripts/figures/figure_S4_PR.pdf', width = 14, height = 10, device = cairo_pdf)
