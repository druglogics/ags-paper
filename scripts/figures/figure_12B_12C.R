library(readr)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(rstatix)

# read file with tumor volume data
tumor_data = readr::read_csv(file = 'data/tumor_vol_data.csv')
tumor_data_wide = tumor_data # keep the wide format for later

# reshape data
tumor_data = tumor_data %>%
  tidyr::pivot_longer(cols = -c(drugs), names_to = 'day', values_to = 'vol') %>%
  mutate(day = as.integer(day)) %>%
  mutate(drugs = factor(x = drugs, levels = c("PI", "Control", "5Z", "5Z-PI")))

pd = position_dodge(0.2)
days = tumor_data %>% distinct(day) %>% pull()

# Figure 12C
tumor_data %>%
  ggpubr::desc_statby(measure.var = 'vol', grps = c('day', 'drugs')) %>%
  ggplot(aes(x = day, y = mean, colour = drugs, group = drugs)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = 'black',
    width = 1, position = pd) +
  geom_line(position = pd) +
  geom_point(position = pd) +
  scale_x_continuous(labels = as.character(days), breaks = days) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'Average tumor volume with SEM', x = 'Days',
    y = latex2exp::TeX('Tumor volume $\\left(mm^3\\right)$')) +
  ylim(c(0, NA)) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())
ggplot2::ggsave(filename = 'scripts/figures/figure_12C.pdf', width = 7, height = 5, device = cairo_pdf)

tumor_diff = tumor_data_wide %>%
  mutate(diff = `19` - `1`, rel_change = (`19`-`1`)/`1`) %>%
  mutate(drugs = factor(x = drugs, levels = c("Control", "PI", "5Z", "5Z-PI"))) %>%
  select(drugs, diff, rel_change)

# Compare single drug vs combo drug group
wilcox_res = tumor_diff %>%
  rstatix::wilcox_test(formula = diff ~ drugs,
    comparisons = list(c('5Z-PI', 'Control'), c('5Z-PI', '5Z'), c('PI','5Z-PI'))) %>%
  select(-`.y.`) %>%
  rstatix::add_xy_position()

# swap heights
y_pos = wilcox_res %>% pull(y.position)
wilcox_res$y.position = y_pos[c(3,1,2)]

# for reproducibility (geom_jitter)
set.seed(42)

# for coloring
set1_col = RColorBrewer::brewer.pal(n = 4, name = 'Set1')

# Figure 12B
tumor_diff %>%
  ggplot(aes(x = drugs, y = diff)) +
  geom_boxplot(aes(fill = drugs)) +
  geom_jitter(position = position_jitter(0.2)) +
  ggpubr::stat_pvalue_manual(wilcox_res, label = "p = {p.adj} ({p.adj.signif})") +
  scale_fill_manual(values = set1_col[c(2,1,3,4)]) +
  labs(title = 'Relative tumor size (Day 1 vs Day 19)', x = "",
    y = latex2exp::TeX('Difference in tumor volume $\\left(mm^3\\right)$')) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())
ggplot2::ggsave(filename = 'scripts/figures/figure_12B.pdf', width = 7, height = 5, device = cairo_pdf)
