library(rstatix)
library(ggpubr)
library(ggplot2)

# load the data
# see script 'scripts/erk_perf_tidy_data.R'
res = readRDS(file = "data/res_erk.rds")

set1_col = RColorBrewer::brewer.pal(9, 'Set1')

# for reproducibility for `add = "jitter"`
set.seed(42)

stat_test_roc = res %>%
  rstatix::wilcox_test(formula = roc_auc ~ erk_state) %>%
  rstatix::add_significance("p")

# ROC AUCs
erk_roc = ggpubr::ggboxplot(res, x = "erk_state", y = "roc_auc", fill = "erk_state",
  palette = c(set1_col[3], set1_col[1]),
  add = "jitter", xlab = "", ylab = "ROC AUC",
  title = "ERK_f Pools Performance Comparison (ROC)") +
  scale_x_discrete(breaks = c("active","inhibited"), labels = c("Active", "Inhibited")) +
  ggpubr::stat_pvalue_manual(stat_test_roc, label = "p = {p} ({p.signif})", y.position = c(1)) +
  ylim(c(0,1)) +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 14)) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 2.1, y = 0.55, label = "Random Predictions (AUC = 0.5)")) +
  theme(legend.position = "none")
ggplot2::ggsave(filename = 'scripts/figures/figure_8_ROC.pdf', width = 7, height = 5, device = cairo_pdf)

stat_test_pr = res %>%
  rstatix::wilcox_test(formula = pr_auc ~ erk_state) %>%
  rstatix::add_significance("p") %>%
  rstatix::add_y_position()

# PR AUCs
erk_pr = ggpubr::ggboxplot(res, x = "erk_state", y = "pr_auc", fill = "erk_state",
  palette = c(set1_col[3], set1_col[1]),
  add = "jitter", xlab = "", ylab = "PR AUC",
  title = "ERK_f Pools Performance Comparison (Precision-Recall)") +
  scale_x_discrete(breaks = c("active","inhibited"), labels = c("Active", "Inhibited")) +
  ggpubr::stat_pvalue_manual(stat_test_pr, label = "p = {p} ({p.signif})") +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 14)) +
  # 6 out 153 drug combos were synergistic
  geom_hline(yintercept = 6/153, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 1.05, y = 0.06, label = "Random Predictions (AUC = 0.04)")) +
  theme(legend.position = "none")
ggplot2::ggsave(filename = 'scripts/figures/figure_8_PR.pdf', width = 7, height = 5, device = cairo_pdf)

# Both plots in one (Figure 8 in publication)
ggpubr::ggarrange(erk_roc, erk_pr, labels = LETTERS[1:2], nrow = 1, ncol = 2)
ggplot2::ggsave(filename = 'scripts/figures/figure_8.pdf', width = 14, height = 5, device = cairo_pdf)
