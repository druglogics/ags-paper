library(ggpubr)
library(ggplot2)
library(scales)

# see `scripts/fit_vs_perf_cascade2_lo.R`
res = readRDS(file = "data/res_fit_aucs.rds")

# data is not normally distributed
shapiro.test(x = res$roc_auc)
shapiro.test(x = res$pr_auc)
shapiro.test(x = res$avg_fit)

# Figure 4 - ROC
p = ggpubr::ggscatter(data = res, x = "avg_fit", y = "roc_auc", color = "avg_fit",
  xlab = "Average Fitness per Model Ensemble",
  title = "Fitness to AGS Steady State vs Performance (ROC)",
  ylab = "ROC AUC", add = "reg.line", conf.int = TRUE,
  add.params = list(color = "blue", fill = "lightgray"),
  cor.coef = TRUE, cor.coeff.args = list(method = "kendall",
    label.y.npc = "bottom", size = 6, cor.coef.name = "tau")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_distiller(n.breaks = 5, labels = scales::label_percent(accuracy = 1L), palette = 'RdYlBu', guide = guide_colourbar(title = '%Data Flipped'))
ggpubr::ggpar(p, legend = "right", font.legend = 14)
ggplot2::ggsave(filename = 'scripts/figures/figure_4_ROC.pdf', width = 7, height = 5, device = cairo_pdf)

# Figure 4 - PR
p = ggpubr::ggscatter(data = res, x = "avg_fit", y = "pr_auc", color = "avg_fit",
  xlab = "Average Fitness per Model Ensemble",
  title = "Fitness to AGS Steady State vs Performance (Precision-Recall)",
  add.params = list(color = "blue", fill = "lightgray"),
  ylab = "PR AUC", add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.coeff.args = list(method = "kendall",
    size = 6, cor.coef.name = "tau")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_distiller(n.breaks = 5, labels = scales::label_percent(accuracy = 1L), palette = 'RdYlBu', guide = guide_colourbar(title = '%Data Flipped'))
ggpubr::ggpar(p, legend = "right", font.legend = 14)
ggplot2::ggsave(filename = 'scripts/figures/figure_4_PR.pdf', width = 7, height = 5, device = cairo_pdf)
