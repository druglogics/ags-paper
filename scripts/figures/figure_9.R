library(rstatix)
library(ggpubr)
library(ggplot2)

# load the data
# see 'scripts/get_param_comp_boot_data.R'
res = readRDS(file = "data/res_param_boot_aucs.rds")

# filter data (compare link mutations vs topology mutation results)
res = res %>%
  filter(param != "topo-and-link")

param_comp = list(c("link-only","topology-only"))

stat_test_roc = res %>%
  rstatix::wilcox_test(formula = roc_auc ~ param, comparisons = param_comp) %>%
  rstatix::add_significance("p")

# ROC AUCs
roc_fig = ggpubr::ggboxplot(res, x = "param", y = "roc_auc", fill = "param", palette = "Set1",
  add = "jitter", xlab = "", ylab = "ROC AUC") +
  scale_x_discrete(breaks = c("link-only","topology-only"),
    labels = c("Parameterization Mutations", "Topology Mutations")) +
  ggpubr::stat_pvalue_manual(stat_test_roc, label = "p = {p} ({p.signif})", y.position = c(1)) +
  ylim(c(0.2,1)) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 16)) +
  geom_hline(yintercept = 0.5, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 2.1, y = 0.45, label="Random Predictions (AUC = 0.5)", size = 5)) +
  theme(legend.position = "none")
ggplot2::ggsave(filename = 'scripts/figures/figure_9_ROC.pdf', width = 7, height = 5, device = cairo_pdf)

stat_test_pr = res %>%
  rstatix::wilcox_test(formula = pr_auc ~ param, comparisons = param_comp) %>%
  rstatix::add_significance("p") %>%
  rstatix::add_y_position()

# PR AUCs (Figure 9)
pr_fig = ggboxplot(res, x = "param", y = "pr_auc", fill = "param", palette = "Set1",
  add = "jitter", xlab = "", ylab = "Precision-Recall AUC") +
  scale_x_discrete(breaks = c("link-only","topology-only"),
    labels = c("Parameterization Mutations", "Topology Mutations")) +
  ggpubr::stat_pvalue_manual(stat_test_pr, label = "p = {p} ({p.signif})") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 16)) +
  geom_hline(yintercept = 6/153, linetype = 'dashed', color = "red") +
  geom_text(aes(x = 1.9, y = 0.08, label="Random Predictions (AUC = 0.04)"), size = 5) +
  theme(legend.position = "none")
ggplot2::ggsave(filename = 'scripts/figures/figure_9_PR.pdf', width = 7, height = 5, device = cairo_pdf)

# Both plots in one
ggpubr::ggarrange(roc_fig, pr_fig, labels = LETTERS[1:2], nrow = 1, ncol = 2)
ggplot2::ggsave(filename = 'scripts/figures/figure_9.pdf', width = 14, height = 5, device = cairo_pdf)
