---
title: "AGS paper I - SI"
author: "[John Zobolas](https://github.com/bblodfon)"
date: "Last updated: 07 May, 2020"
description: "AGS paper I - SI"
url: 'https\://username.github.io/reponame/'
github-repo: "username/reponame"
bibliography: ["references.bib", "packages.bib"]
link-citations: true
site: bookdown::bookdown_site
---



# Intro {-}

This report is the **supplementary material** for the AGS I Paper and has all the simulation results and investigations related to that paper, as well as guidelines for reproducing the results.

A list of things that change between the simulations and the presented graphs are:

- The number of `Gitsbe` simulations: more simulations, more models generated.
- The type of mutation that `Gitsbe` models have:
unless otherwise specified, the `Gitsbe` models have only [link operator mutations](https://druglogics.github.io/druglogics-doc/gitsbe-config.html#genetic-algorithm).
[Topology mutations] were also tested as well as a combination of [topology and link operator mutations].
- The [training data](https://druglogics.github.io/druglogics-doc/training-data.html) for the `Gitsbe` models: *steady state* (calibrated models) vs *proliferative profile* (proliferative models).
Also *randomly generated models* were produced for the link operator mutations using the [abmlog module](https://github.com/druglogics/abmlog).
- The type of mathematical model (HSA or Bliss) used in `Drabme` to evaluate the synergies either from the [@Flobak2015] for the [Cascade 1.0 Analysis] or from the [@Flobak2019] dataset for the [Cascade 2.0 Analysis].
More info on the calcualtions that Drabme does [see here](https://druglogics.github.io/druglogics-doc/drabme-description.html#drabme-description).
- The type of output used from `Drabme`: ensemble-wise or model-wise [synergy results](https://druglogics.github.io/druglogics-doc/drabme-install.html#drabme-output).

For the ROC curves we used the function `get_roc_stats()` from [@R-usefun] and for the PR curves the `pr.curve()` from [@R-PRROC] (see also [@Grau2015]).

The report template is from the `rtemps` R package [@R-rtemps].

# Input {-}

Loading libraries:

```r
library(DT)
library(ggpubr)
library(RColorBrewer)
library(xfun)
library(dplyr)
library(tibble)
library(emba)
library(usefun)
library(readr)
library(stringr)
library(latex2exp)
library(corrplot)
library(PRROC)
library(equatiomatic)
library(glmnet)
library(knitr)
library(MAMSE)
```

# Cascade 1.0 Analysis {-}

:::{.blue-box}
Performance of automatically parameterized models against published data in [@Flobak2015]
:::

## HSA results {-}

:::{.note}
- *HSA* refers to the synergy method used in `Drabme` to assess the synergies from the `gitsbe` models
- We test performance using ROC and PR AUC for both the *ensemble-wise* and *model-wise* synergies from `Drabme`
- **Calibrated** models: fitted to steady state ($50$ simulations)
- **Proliferative** models: fitted to [proliferation profile](https://druglogics.github.io/druglogics-doc/training-data.html#unperturbed-condition---globaloutput-response)
- **Random** models: produced via `abmlog` (see [here](#random-model-results) and used in `Drabme` with `synergy_method: hsa`
- Mutations on **link operator** only
:::

Load results:

```r
# 'ss' => calibrated models, 'random' => random models, 'prolif' => proliferative models
# 'ew' => ensemble-wise, 'mw' => model-wise

## HSA results
ss_hsa_ew_file = paste0("results/hsa/cascade_1.0_ss_50sim_fixpoints_ensemblewise_synergies.tab")
ss_hsa_mw_file = paste0("results/hsa/cascade_1.0_ss_50sim_fixpoints_modelwise_synergies.tab")
prolif_hsa_ew_file = paste0("results/hsa/cascade_1.0_rand_50sim_fixpoints_ensemblewise_synergies.tab")
prolif_hsa_mw_file = paste0("results/hsa/cascade_1.0_rand_50sim_fixpoints_modelwise_synergies.tab")
random_hsa_ew_file = paste0("results/hsa/cascade_1.0_random_ensemblewise_synergies.tab")
random_hsa_mw_file = paste0("results/hsa/cascade_1.0_random_modelwise_synergies.tab")

ss_hsa_ensemblewise_synergies = emba::get_synergy_scores(ss_hsa_ew_file)
ss_hsa_modelwise_synergies = emba::get_synergy_scores(ss_hsa_mw_file, file_type = "modelwise")
prolif_hsa_ensemblewise_synergies = emba::get_synergy_scores(prolif_hsa_ew_file)
prolif_hsa_modelwise_synergies = emba::get_synergy_scores(prolif_hsa_mw_file, file_type = "modelwise")
random_hsa_ensemblewise_synergies = emba::get_synergy_scores(random_hsa_ew_file)
random_hsa_modelwise_synergies = emba::get_synergy_scores(random_hsa_mw_file, file_type = "modelwise")

# calculate probability of synergy in the modelwise results
ss_hsa_modelwise_synergies = ss_hsa_modelwise_synergies %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
prolif_hsa_modelwise_synergies = prolif_hsa_modelwise_synergies %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
random_hsa_modelwise_synergies = random_hsa_modelwise_synergies %>%
  mutate(synergy_prob_random = synergies/(synergies + `non-synergies`))

observed_synergies_file = paste0("results/observed_synergies_cascade_1.0")
observed_synergies = get_observed_synergies(observed_synergies_file)
# 1 (positive/observed synergy) or 0 (negative/not observed) for all tested drug combinations
observed = sapply(random_hsa_modelwise_synergies$perturbation %in% observed_synergies, as.integer)
```

### ROC curves {-}


```r
# 'ew' => ensemble-wise, 'mw' => model-wise
pred_ew_hsa = bind_cols(ss_hsa_ensemblewise_synergies %>% rename(ss_score = score),
  prolif_hsa_ensemblewise_synergies %>% select(score) %>% rename(prolif_score = score),
  random_hsa_ensemblewise_synergies %>% select(score) %>% rename(random_score = score), 
  as_tibble_col(observed, column_name = "observed"))

pred_mw_hsa = bind_cols(
  ss_hsa_modelwise_synergies %>% select(perturbation, synergy_prob_ss),
  prolif_hsa_modelwise_synergies %>% select(synergy_prob_ss) %>% rename(synergy_prob_prolif = synergy_prob_ss),
  random_hsa_modelwise_synergies %>% select(synergy_prob_random),
  as_tibble_col(observed, column_name = "observed"))

res_ss_ew = get_roc_stats(df = pred_ew_hsa, pred_col = "ss_score", label_col = "observed")
res_prolif_ew = get_roc_stats(df = pred_ew_hsa, pred_col = "prolif_score", label_col = "observed")
res_random_ew = get_roc_stats(df = pred_ew_hsa, pred_col = "random_score", label_col = "observed")

res_ss_mw = get_roc_stats(df = pred_mw_hsa, pred_col = "synergy_prob_ss", label_col = "observed", direction = ">")
res_prolif_mw = get_roc_stats(df = pred_mw_hsa, pred_col = "synergy_prob_prolif", label_col = "observed", direction = ">")
res_random_mw = get_roc_stats(df = pred_mw_hsa, pred_col = "synergy_prob_random", label_col = "observed", direction = ">")

# Plot ROCs
my_palette = RColorBrewer::brewer.pal(n = 9, name = "Set1")

plot(x = res_ss_ew$roc_stats$FPR, y = res_ss_ew$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (HSA)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_prolif_ew$roc_stats$FPR, y = res_prolif_ew$roc_stats$TPR, 
  lwd = 3, col = my_palette[2])
lines(x = res_random_ew$roc_stats$FPR, y = res_random_ew$roc_stats$TPR, 
  lwd = 3, col = my_palette[3])
legend('bottomright', title = 'AUC', col = my_palette[1:3], pch = 19,
  legend = c(paste(round(res_ss_ew$AUC, digits = 3), "Calibrated"), 
    paste(round(res_prolif_ew$AUC, digits = 3), "Proliferative"),
    paste(round(res_random_ew$AUC, digits = 3), "Random")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(x = res_ss_mw$roc_stats$FPR, y = res_ss_mw$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Model-wise synergies (HSA)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_prolif_mw$roc_stats$FPR, y = res_prolif_mw$roc_stats$TPR, 
  lwd = 3, col = my_palette[2])
lines(x = res_random_mw$roc_stats$FPR, y = res_random_mw$roc_stats$TPR, 
  lwd = 3, col = my_palette[3])
legend('bottomright', title = 'AUC', col = my_palette[1:3], pch = 19,
  legend = c(paste(round(res_ss_mw$AUC, digits = 3), "Calibrated"),
    paste(round(res_prolif_mw$AUC, digits = 3), "Proliferative"), 
    paste(round(res_random_mw$AUC, digits = 3), "Random")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
```

<img src="index_files/figure-html/ROC HSA Cascade 1.0-1.png" width="50%" /><img src="index_files/figure-html/ROC HSA Cascade 1.0-2.png" width="50%" />

### PR curves {-}


```r
pr_ss_ew_hsa = pr.curve(scores.class0 = pred_ew_hsa %>% pull(ss_score) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_prolif_ew_hsa = pr.curve(scores.class0 = pred_ew_hsa %>% pull(prolif_score) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_hsa %>% pull(observed), curve = TRUE)
pr_random_ew_hsa = pr.curve(scores.class0 = pred_ew_hsa %>% pull(random_score) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_hsa %>% pull(observed), curve = TRUE)

pr_ss_mw_hsa = pr.curve(scores.class0 = pred_mw_hsa %>% pull(synergy_prob_ss), 
  weights.class0 = pred_mw_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_prolif_mw_hsa = pr.curve(scores.class0 = pred_mw_hsa %>% pull(synergy_prob_prolif), 
  weights.class0 = pred_mw_hsa %>% pull(observed), curve = TRUE)
pr_random_mw_hsa = pr.curve(scores.class0 = pred_mw_hsa %>% pull(synergy_prob_random), 
  weights.class0 = pred_mw_hsa %>% pull(observed), curve = TRUE)

plot(pr_ss_ew_hsa, main = 'PR curve, Ensemble-wise synergies (HSA)', 
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_prolif_ew_hsa, add = TRUE, color = my_palette[2])
plot(pr_random_ew_hsa, add = TRUE, color = my_palette[3])
legend('topright', title = 'AUC', col = my_palette[1:3], pch = 19,
  legend = c(paste(round(pr_ss_ew_hsa$auc.davis.goadrich, digits = 3), "Calibrated"), 
    paste(round(pr_prolif_ew_hsa$auc.davis.goadrich, digits = 3), "Proliferative"),
    paste(round(pr_random_ew_hsa$auc.davis.goadrich, digits = 3), "Random")))
grid(lwd = 0.5)

plot(pr_ss_mw_hsa, main = 'PR curve, Model-wise synergies (HSA)', 
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_prolif_mw_hsa, add = TRUE, color = my_palette[2])
plot(pr_random_mw_hsa, add = TRUE, color = my_palette[3])
legend('left', title = 'AUC', col = my_palette[1:3], pch = 19,
  legend = c(paste(round(pr_ss_mw_hsa$auc.davis.goadrich, digits = 3), "Calibrated"),
    paste(round(pr_prolif_mw_hsa$auc.davis.goadrich, digits = 3), "Proliferative"),
    paste(round(pr_random_mw_hsa$auc.davis.goadrich, digits = 3), "Random")))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/PR HSA Cascade 1.0-1.png" width="50%" /><img src="index_files/figure-html/PR HSA Cascade 1.0-2.png" width="50%" />

:::{.green-box}
- Calibrated models perform a lot better than either the random or proliferative ones
:::

### AUC sensitivity {-}

:::{#auc-sensitivity .blue-box}
- Investigate **combining the synergy results of calibrated and proliferative (random) models**
- Quantify the amount of information from the *proliferative* (*random*) models that can be used to augment the calibrated results?
- **Ensemble-wise** scenario: $score = calibrated + \beta \times proliferative$ ($random$)
  - $\beta \rightarrow +\infty$: mostly **proliferative (random) model predictions**
  - $\beta \rightarrow -\infty$: mostly **reverse proliferative (random) model predictions**
- **Model-wise** scenario: $(1-w) \times prob_{cal} + w \times prob_{rand}, w \in[0,1]$
  - $w=0$: only calibrated model predictions
  - $w=1$: only proliferative (random) model predictions
:::


```r
# Ensemble-wise
betas = seq(from = -20, to = 20, by = 0.1)

random_roc = sapply(betas, function(beta) {
  pred_ew_hsa = pred_ew_hsa %>% mutate(combined_score = ss_score + beta * random_score)
  res = roc.curve(scores.class0 = pred_ew_hsa %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_ew_hsa %>% pull(observed))
  auc_value = res$auc
})

prolif_roc = sapply(betas, function(beta) {
  pred_ew_hsa = pred_ew_hsa %>% mutate(combined_score = ss_score + beta * prolif_score)
  res = roc.curve(scores.class0 = pred_ew_hsa %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_ew_hsa %>% pull(observed))
  auc_value = res$auc
})

random_pr = sapply(betas, function(beta) {
  pred_ew_hsa = pred_ew_hsa %>% mutate(combined_score = ss_score + beta * random_score)
  res = pr.curve(scores.class0 = pred_ew_hsa %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_ew_hsa %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

prolif_pr = sapply(betas, function(beta) {
  pred_ew_hsa = pred_ew_hsa %>% mutate(combined_score = ss_score + beta * prolif_score)
  res = pr.curve(scores.class0 = pred_ew_hsa %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_ew_hsa %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_ew = as_tibble(cbind(betas, random_roc, prolif_roc, random_pr, prolif_pr))
df_ew = df_ew %>% tidyr::pivot_longer(-betas, names_to = "type", values_to = "AUC")

ggline(data = df_ew, x = "betas", y = "AUC", numeric.x.axis = TRUE, color = "type",
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under Curve)", 
  legend = "none", facet.by = "type", palette = my_palette,
  panel.labs = list(type = c("PR: calibrated + β x proliferative", 
    "ROC: calibrated + β x proliferative", "PR: calibrated + β x random", 
    "ROC: calibrated + β x random")),
  title = TeX("AUC sensitivity to $\\beta$ parameter (HSA, Cascade 1.0)")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = -1, color = "red", size = 0.3, linetype = "dashed") + 
  geom_text(aes(x=-2, label="β = -1", y=0.25), colour="black", angle=90) + 
  grids()
```

<img src="index_files/figure-html/AUC sensitivity (HSA - Cascade 1.0)-1.png" width="2100" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

random_roc_mw = sapply(weights, function(w) {
  pred_mw_hsa = pred_mw_hsa %>%
    mutate(weighted_prob = (1 - w) * pred_mw_hsa$synergy_prob_ss + w * pred_mw_hsa$synergy_prob_random)
  res = roc.curve(scores.class0 = pred_mw_hsa %>% pull(weighted_prob),
    weights.class0 = pred_mw_hsa %>% pull(observed))
  auc_value = res$auc
})

prolif_roc_mw = sapply(weights, function(w) {
  pred_mw_hsa = pred_mw_hsa %>%
    mutate(weighted_prob = (1 - w) * pred_mw_hsa$synergy_prob_ss + w * pred_mw_hsa$synergy_prob_prolif)
  res = roc.curve(scores.class0 = pred_mw_hsa %>% pull(weighted_prob),
    weights.class0 = pred_mw_hsa %>% pull(observed))
  auc_value = res$auc
})

random_pr_mw = sapply(weights, function(w) {
  pred_mw_hsa = pred_mw_hsa %>% 
    mutate(weighted_prob = (1 - w) * pred_mw_hsa$synergy_prob_ss + w * pred_mw_hsa$synergy_prob_random)
  res = pr.curve(scores.class0 = pred_mw_hsa %>% pull(weighted_prob), 
    weights.class0 = pred_mw_hsa %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

prolif_pr_mw = sapply(weights, function(w) {
  pred_mw_hsa = pred_mw_hsa %>% 
    mutate(weighted_prob = (1 - w) * pred_mw_hsa$synergy_prob_ss + w * pred_mw_hsa$synergy_prob_prolif)
  res = pr.curve(scores.class0 = pred_mw_hsa %>% pull(weighted_prob), 
    weights.class0 = pred_mw_hsa %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_mw = as_tibble(cbind(weights, random_roc_mw, prolif_roc_mw, random_pr_mw, prolif_pr_mw))
df_mw = df_mw %>% tidyr::pivot_longer(-weights, names_to = "type", values_to = "AUC")

ggline(data = df_mw, x = "weights", y = "AUC", numeric.x.axis = TRUE, color = "type",
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "AUC (Area Under Curve)", 
  legend = "none", facet.by = "type", palette = my_palette,
  panel.labs = list(type = c("PR: (1-w) x prob(ss) + w x prob(prolif)", 
    "ROC: (1-w) x prob(ss) + w x prob(prolif)", "PR: (1-w) x prob(ss) + w x prob(random)", 
    "ROC: (1-w) x prob(ss) + w x prob(random)")), title.position = "center",
  title = TeX("AUC sensitivity to weighted average score (HSA, Cascade 1.0)")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  grids()
```

<img src="index_files/figure-html/AUC sensitivity (HSA - Cascade 1.0)-2.png" width="2100" style="display: block; margin: auto;" />

:::{.green-box}
- There are $\beta$ values that can boost the predictive performance of the calibrated models (ensemble-wise) but no $w$ weight in the model-wise case.
- $\beta=-1$ seems to be a common value that maximizes both the ROC-AUC and the PR-AUC.
- The PR-AUC is **more sensitive** than the ROC-AUC, so a better indicator of performance.
:::

## Bliss results {-}

:::{.note}
- *Bliss* refers to the synergy method used in `Drabme` to assess the synergies from the `gitsbe` models
- We test performance using ROC AUC for both the *ensemble-wise* and *model-wise* synergies from `Drabme`
- **Calibrated** models: fitted to steady state
- **Random** models: produced via `abmlog` (see [here](#random-model-results) and used in `Drabme` with `synergy_method: bliss`
:::

Load results:

```r
# 'ss' => calibrated models, 'random' => random models

## Bliss results
ss_bliss_ensemblewise_file = paste0("results/bliss/cascade_1.0_ss_50sim_fixpoints_ensemblewise_synergies.tab")
ss_bliss_modelwise_file = paste0("results/bliss/cascade_1.0_ss_50sim_fixpoints_modelwise_synergies.tab")
random_bliss_ensemblewise_file = paste0("results/bliss/cascade_1.0_random_bliss_ensemblewise_synergies.tab")
random_bliss_modelwise_file = paste0("results/bliss/cascade_1.0_random_bliss_modelwise_synergies.tab")

ss_bliss_ensemblewise_synergies = emba::get_synergy_scores(ss_bliss_ensemblewise_file)
ss_bliss_modelwise_synergies = emba::get_synergy_scores(ss_bliss_modelwise_file, file_type = "modelwise")
random_bliss_ensemblewise_synergies = emba::get_synergy_scores(random_bliss_ensemblewise_file)
random_bliss_modelwise_synergies = emba::get_synergy_scores(random_bliss_modelwise_file, file_type = "modelwise")

# calculate probability of synergy in the modelwise results
ss_bliss_modelwise_synergies = ss_bliss_modelwise_synergies %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
random_bliss_modelwise_synergies = random_bliss_modelwise_synergies %>%
  mutate(synergy_prob_random = synergies/(synergies + `non-synergies`))
```

### ROC curves {-}


```r
# 'ew' => ensemble-wise, 'mw' => model-wise
pred_ew_bliss = bind_cols(ss_bliss_ensemblewise_synergies %>% rename(ss_score = score), 
  random_bliss_ensemblewise_synergies %>% select(score) %>% rename(random_score = score), 
  as_tibble_col(observed, column_name = "observed"))

pred_mw_bliss = bind_cols(
  ss_bliss_modelwise_synergies %>% select(perturbation, synergy_prob_ss),
  random_bliss_modelwise_synergies %>% select(synergy_prob_random),
  as_tibble_col(observed, column_name = "observed"))

res_ss_ew = get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score", label_col = "observed")
res_random_ew = get_roc_stats(df = pred_ew_bliss, pred_col = "random_score", label_col = "observed")

res_ss_mw = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_ss", label_col = "observed", direction = ">")
res_random_mw = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_random", label_col = "observed", direction = ">")

# Plot ROCs
plot(x = res_ss_ew$roc_stats$FPR, y = res_ss_ew$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_random_ew$roc_stats$FPR, y = res_random_ew$roc_stats$TPR, 
  lwd = 3, col = my_palette[2])
legend('bottomright', title = 'AUC', col = my_palette[1:2], pch = 19,
  legend = c(paste(round(res_ss_ew$AUC, digits = 3), "Calibrated"), 
    paste(round(res_random_ew$AUC, digits = 3), "Random")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(x = res_ss_mw$roc_stats$FPR, y = res_ss_mw$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Model-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_random_mw$roc_stats$FPR, y = res_random_mw$roc_stats$TPR, 
  lwd = 3, col = my_palette[2])
legend('bottomright', title = 'AUC', col = my_palette[1:2], pch = 19,
  legend = c(paste(round(res_ss_mw$AUC, digits = 3), "Calibrated"), 
    paste(round(res_random_mw$AUC, digits = 3), "Random")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
```

<img src="index_files/figure-html/ROC Bliss Cascade 1.0-1.png" width="50%" /><img src="index_files/figure-html/ROC Bliss Cascade 1.0-2.png" width="50%" />

### PR curves {-}


```r
pr_ss_ew_bliss = pr.curve(scores.class0 = pred_ew_bliss %>% pull(ss_score) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_random_ew_bliss = pr.curve(scores.class0 = pred_ew_bliss %>% pull(random_score) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)

pr_ss_mw_bliss = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_ss), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_random_mw_bliss = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_random), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE)

plot(pr_ss_ew_bliss, main = 'PR curve, Ensemble-wise synergies (Bliss)', 
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_random_ew_bliss, add = TRUE, color = my_palette[2])
legend('bottomright', title = 'AUC', col = my_palette[1:2], pch = 19, cex = 0.8,
  legend = c(paste(round(pr_ss_ew_bliss$auc.davis.goadrich, digits = 3), "Calibrated"), 
    paste(round(pr_random_ew_bliss$auc.davis.goadrich, digits = 3), "Random")))
grid(lwd = 0.5)

plot(pr_ss_mw_bliss, main = 'PR curve, Model-wise synergies (Bliss)', 
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_random_mw_bliss, add = TRUE, color = my_palette[2])
legend('left', title = 'AUC', col = my_palette[1:2], pch = 19,
  legend = c(paste(round(pr_ss_mw_bliss$auc.davis.goadrich, digits = 3), "Calibrated"), 
    paste(round(pr_random_mw_bliss$auc.davis.goadrich, digits = 3), "Random")))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/PR Bliss Cascade 1.0-1.png" width="50%" /><img src="index_files/figure-html/PR Bliss Cascade 1.0-2.png" width="50%" />

### ROC AUC sensitivity {-}

Investigate same thing as described in [here](#auc-sensitivity).


```r
# Ensemble-wise
betas = seq(from = -20, to = 20, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_ew_bliss = pred_ew_bliss %>% mutate(combined_score = ss_score + beta * random_score)
  res = get_roc_stats(df = pred_ew_bliss, pred_col = "combined_score", label_col = "observed")
  auc_value = res$AUC
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2]) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/AUC sensitivity (Bliss - cascade 1.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

auc_values_mw = sapply(weights, function(w) {
  pred_mw_bliss = pred_mw_bliss %>% 
    mutate(weighted_prob = (1 - w) * pred_mw_bliss$synergy_prob_ss + w * pred_mw_bliss$synergy_prob_random)
  res = get_roc_stats(df = pred_mw_bliss, pred_col = "weighted_prob", label_col = "observed", direction = ">")
  auc_value = res$AUC
})

df_mw = as_tibble(cbind(weights, auc_values_mw))

ggline(data = df_mw, x = "weights", y = "auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/AUC sensitivity (Bliss - cascade 1.0)-2.png" width="80%" style="display: block; margin: auto;" />

:::{.green-box}
- Symmetricity (Ensemble-wise): $AUC_{\beta \rightarrow +\infty} + AUC_{\beta \rightarrow -\infty} \approx 1$
- Random models perform worse than calibrated ones
- There are $\beta$ values that can boost the predictive performance of the combined synergy classifier but no $w$ weight in the model-wise case
:::

### PR AUC sensitivity {-}


```r
# Ensemble-wise
betas = seq(from = -20, to = 20, by = 0.1)

pr_auc_values_ew = sapply(betas, function(beta) {
  pred_ew_bliss = pred_ew_bliss %>% mutate(combined_score = ss_score + beta * random_score)
  res = pr.curve(scores.class0 = pred_ew_bliss %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_ew_bliss %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_pr_ew = as_tibble(cbind(betas, pr_auc_values_ew))

ggline(data = df_pr_ew, x = "betas", y = "pr_auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "PR-AUC (Area Under PR Curve)",
  title = TeX("PR-AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2]) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/PR AUC sensitivity (Bliss - Cascade 1.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

pr_auc_values_mw = sapply(weights, function(w) {
  pred_mw_bliss = pred_mw_bliss %>% 
    mutate(weighted_prob = (1 - w) * pred_mw_bliss$synergy_prob_ss + w * pred_mw_bliss$synergy_prob_random)
  res = pr.curve(scores.class0 = pred_mw_bliss %>% pull(weighted_prob), 
    weights.class0 = pred_mw_bliss %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_pr_mw = as_tibble(cbind(weights, pr_auc_values_mw))

ggline(data = df_pr_mw, x = "weights", y = "pr_auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "PR-AUC (Area Under PR Curve)",
  title = TeX("PR-AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/PR AUC sensitivity (Bliss - Cascade 1.0)-2.png" width="80%" style="display: block; margin: auto;" />

## Correlation {-}

We test for correlation between all the results shown in the ROC curves.
This means *ensemble-wise* vs *model-wise*, *random* models vs *calibrated (ss)* models and *HSA* vs *Bliss* synergy assessment.
*P-values* are represented at 3 significant levels: $0.05, 0.01, 0.001$ (\*, \*\*, \*\*\*) and the correlation coefficient is calculated using Kendall's *tau* statistic.


```r
synergy_scores = bind_cols(
  pred_ew_hsa %>% select(ss_score, random_score) %>% rename(ss_ensemble_hsa = ss_score, random_ensemble_hsa = random_score),
  pred_ew_bliss %>% select(ss_score, random_score) %>% rename(ss_ensemble_bliss = ss_score, random_ensemble_bliss = random_score),
  pred_mw_hsa %>% select(synergy_prob_ss, synergy_prob_random) %>% 
    rename(ss_modelwise_hsa = synergy_prob_ss, random_modelwise_hsa = synergy_prob_random),
  pred_mw_bliss %>% select(synergy_prob_ss, synergy_prob_random) %>% 
    rename(ss_modelwise_bliss = synergy_prob_ss, random_modelwise_bliss = synergy_prob_random)
  )

M = cor(synergy_scores, method = "kendall")
res = cor.mtest(synergy_scores, method = "kendall")
corrplot(corr = M, type = "upper", p.mat = res$p, sig.level = c(.001, .01, .05), 
  pch.cex = 1, pch.col = "white", insig = "label_sig", tl.col = "black", tl.srt = 45)
```

<img src="index_files/figure-html/Correlation of ROC results (Cascade 1.0)-1.png" width="2100" />

:::{.green-box}
- **HSA and Bliss results correlate**, higher for the model-wise than the ensemble-wise results.
- **Model-wise don't correlate with ensemble-wise results**.
:::

## Fitness Evolution {-}

Results are from the simulation result with $50$ Gitsbe simulations, fitting to steady state (**calibrated models**) and *HSA* Drabme synergy assessment.
We show only $10$ simulations - the first ones that spanned the maximum defined generations in the configuration ($20$), meaning that they did not surpass the target fitness threhold specified ($0.99$).
Each data point is the average fitness in that generation out of $20$ models.


```r
fitness_summary_file = paste0("results/hsa/cascade_1.0_ss_50sim_fixpoints_summary.txt")

read_summary_file = function(file_name) {
  lines = readr::read_lines(file = fitness_summary_file, skip = 5, skip_empty_rows = TRUE)
  
  data_list = list()
  index = 1
  
  gen_fit_list = list()
  gen_index = 1
  for (line_index in 1:length(lines)) {
    line = lines[line_index]
    if (stringr::str_detect(string = line, pattern = "Simulation")) {
      data_list[[index]] = bind_cols(gen_fit_list)
      index = index + 1
      
      gen_fit_list = list()
      gen_index = 1
    } else { # read fitness values
      gen_fit_list[[gen_index]] = as_tibble_col(as.numeric(unlist(strsplit(line, split = '\t'))))
      gen_index = gen_index + 1
    }
  }
  
  # add the last simulation's values
  data_list[[index]] = bind_cols(gen_fit_list)
  
  return(data_list)
}

fit_res = read_summary_file(file_name = fitness_summary_file)

first_sim_data = colMeans(fit_res[[1]])
plot(1:length(first_sim_data), y = first_sim_data, ylim = c(0,1), 
  xlim = c(0,20), type = 'l', lwd = 1.5, 
  main = 'Fitness vs Generation (10 Simulations)', xlab = 'Generations', 
  ylab = 'Average Fitness', col = usefun:::colors.100[1])
index = 2
for (fit_data in fit_res) {
  if (index > 10) break
  if (ncol(fit_data) != 20) next
  mean_fit_per_gen = colMeans(fit_data)
  lines(x = 1:length(mean_fit_per_gen), y = mean_fit_per_gen, lwd = 1.5,
    col = usefun:::colors.100[index])
  index = index + 1
}
grid(lwd = 0.5)
```

<img src="index_files/figure-html/fitness evolution-1.png" width="2100" style="display: block; margin: auto;" />


# Cascade 2.0 Analysis {-}

:::{.blue-box}
Performance of automatically parameterized models against a new dataset (SINTEF, AGS only)
:::

## Calibrated vs Random (HSA) {-}

:::{.note}
- *HSA* refers to the synergy method used in `Drabme` to assess the synergies from the `gitsbe` models
- We test performance using ROC AUC for both the *ensemble-wise* and *model-wise* synergies from `Drabme`
- **Calibrated** models: fitted to steady state
- **Random** models: produced via `abmlog` (see [here](#random-model-results) and used in `Drabme` with `synergy_method: hsa`
:::

Load results:


```r
# 'ss' => calibrated models, 'random' => random models

## HSA results
ss_hsa_ensemblewise_50sim_file = paste0("results/hsa/cascade_2.0_ss_50sim_fixpoints_ensemblewise_synergies.tab")
ss_hsa_modelwise_50sim_file = paste0("results/hsa/cascade_2.0_ss_50sim_fixpoints_modelwise_synergies.tab")
ss_hsa_ensemblewise_100sim_file = paste0("results/hsa/cascade_2.0_ss_100sim_fixpoints_ensemblewise_synergies.tab")
ss_hsa_modelwise_100sim_file = paste0("results/hsa/cascade_2.0_ss_100sim_fixpoints_modelwise_synergies.tab")
ss_hsa_ensemblewise_150sim_file = paste0("results/hsa/cascade_2.0_ss_150sim_fixpoints_ensemblewise_synergies.tab")
ss_hsa_modelwise_150sim_file = paste0("results/hsa/cascade_2.0_ss_150sim_fixpoints_modelwise_synergies.tab")
ss_hsa_ensemblewise_200sim_file = paste0("results/hsa/cascade_2.0_ss_200sim_fixpoints_ensemblewise_synergies.tab")
ss_hsa_modelwise_200sim_file = paste0("results/hsa/cascade_2.0_ss_200sim_fixpoints_modelwise_synergies.tab")
random_hsa_ensemblewise_file = paste0("results/hsa/cascade_2.0_random_ensemblewise_synergies.tab")
random_hsa_modelwise_file = paste0("results/hsa/cascade_2.0_random_modelwise_synergies.tab")

ss_hsa_ensemblewise_synergies_50sim = emba::get_synergy_scores(ss_hsa_ensemblewise_50sim_file)
ss_hsa_modelwise_synergies_50sim = emba::get_synergy_scores(ss_hsa_modelwise_50sim_file, file_type = "modelwise")
ss_hsa_ensemblewise_synergies_100sim = emba::get_synergy_scores(ss_hsa_ensemblewise_100sim_file)
ss_hsa_modelwise_synergies_100sim = emba::get_synergy_scores(ss_hsa_modelwise_100sim_file, file_type = "modelwise")
ss_hsa_ensemblewise_synergies_150sim = emba::get_synergy_scores(ss_hsa_ensemblewise_150sim_file)
ss_hsa_modelwise_synergies_150sim = emba::get_synergy_scores(ss_hsa_modelwise_150sim_file, file_type = "modelwise")
ss_hsa_ensemblewise_synergies_200sim = emba::get_synergy_scores(ss_hsa_ensemblewise_200sim_file)
ss_hsa_modelwise_synergies_200sim = emba::get_synergy_scores(ss_hsa_modelwise_200sim_file, file_type = "modelwise")
random_hsa_ensemblewise_synergies = emba::get_synergy_scores(random_hsa_ensemblewise_file)
random_hsa_modelwise_synergies = emba::get_synergy_scores(random_hsa_modelwise_file, file_type = "modelwise")

# calculate probability of synergy in the modelwise results
ss_hsa_modelwise_synergies_50sim = ss_hsa_modelwise_synergies_50sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
ss_hsa_modelwise_synergies_100sim = ss_hsa_modelwise_synergies_100sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
ss_hsa_modelwise_synergies_150sim = ss_hsa_modelwise_synergies_150sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
ss_hsa_modelwise_synergies_200sim = ss_hsa_modelwise_synergies_200sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
random_hsa_modelwise_synergies = random_hsa_modelwise_synergies %>%
  mutate(synergy_prob_random = synergies/(synergies + `non-synergies`))

observed_synergies_file = paste0("results/observed_synergies_cascade_2.0")
observed_synergies = get_observed_synergies(observed_synergies_file)
# 1 (positive/observed synergy) or 0 (negative/not observed) for all tested drug combinations
observed = sapply(random_hsa_modelwise_synergies$perturbation %in% observed_synergies, as.integer)
```

### ROC curves {-}


```r
# 'ew' => ensemble-wise, 'mw' => model-wise
pred_ew_hsa = bind_cols(
  ss_hsa_ensemblewise_synergies_50sim %>% select(score) %>% rename(ss_score_50sim = score),
  ss_hsa_ensemblewise_synergies_100sim %>% select(score) %>% rename(ss_score_100sim = score),
  ss_hsa_ensemblewise_synergies_150sim %>% select(score) %>% rename(ss_score_150sim = score),
  ss_hsa_ensemblewise_synergies_200sim %>% select(score) %>% rename(ss_score_200sim = score),
  random_hsa_ensemblewise_synergies %>% select(score) %>% rename(random_score = score), 
  as_tibble_col(observed, column_name = "observed"))

pred_mw_hsa = bind_cols(
  ss_hsa_modelwise_synergies_50sim %>% select(perturbation, synergy_prob_ss) %>% rename(synergy_prob_ss_50sim = synergy_prob_ss),
  ss_hsa_modelwise_synergies_100sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_100sim = synergy_prob_ss),
  ss_hsa_modelwise_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_150sim = synergy_prob_ss),
  ss_hsa_modelwise_synergies_200sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_200sim = synergy_prob_ss),
  random_hsa_modelwise_synergies %>% select(synergy_prob_random),
  as_tibble_col(observed, column_name = "observed"))

res_ss_ew_50sim = get_roc_stats(df = pred_ew_hsa, pred_col = "ss_score_50sim", label_col = "observed")
res_ss_ew_100sim = get_roc_stats(df = pred_ew_hsa, pred_col = "ss_score_100sim", label_col = "observed")
res_ss_ew_150sim = get_roc_stats(df = pred_ew_hsa, pred_col = "ss_score_150sim", label_col = "observed")
res_ss_ew_200sim = get_roc_stats(df = pred_ew_hsa, pred_col = "ss_score_200sim", label_col = "observed")
res_random_ew = get_roc_stats(df = pred_ew_hsa, pred_col = "random_score", label_col = "observed")

res_ss_mw_50sim = get_roc_stats(df = pred_mw_hsa, pred_col = "synergy_prob_ss_50sim", label_col = "observed", direction = ">")
res_ss_mw_100sim = get_roc_stats(df = pred_mw_hsa, pred_col = "synergy_prob_ss_100sim", label_col = "observed", direction = ">")
res_ss_mw_150sim = get_roc_stats(df = pred_mw_hsa, pred_col = "synergy_prob_ss_150sim", label_col = "observed", direction = ">")
res_ss_mw_200sim = get_roc_stats(df = pred_mw_hsa, pred_col = "synergy_prob_ss_200sim", label_col = "observed", direction = ">")
res_random_mw = get_roc_stats(df = pred_mw_hsa, pred_col = "synergy_prob_random", label_col = "observed", direction = ">")

# Plot ROCs
my_palette = RColorBrewer::brewer.pal(n = 9, name = "Set1")

plot(x = res_ss_ew_50sim$roc_stats$FPR, y = res_ss_ew_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (HSA)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_ss_ew_100sim$roc_stats$FPR, y = res_ss_ew_100sim$roc_stats$TPR, 
  lwd = 3, col = my_palette[2])
lines(x = res_ss_ew_150sim$roc_stats$FPR, y = res_ss_ew_150sim$roc_stats$TPR, 
  lwd = 3, col = my_palette[3])
lines(x = res_ss_ew_200sim$roc_stats$FPR, y = res_ss_ew_200sim$roc_stats$TPR, 
  lwd = 3, col = my_palette[4])
lines(x = res_random_ew$roc_stats$FPR, y = res_random_ew$roc_stats$TPR, 
  lwd = 3, col = my_palette[5])
legend('bottomright', title = 'AUC', col = my_palette[1:5], pch = 19,
  legend = c(paste(round(res_ss_ew_50sim$AUC, digits = 3), "Calibrated (50 sim)"),
    paste(round(res_ss_ew_100sim$AUC, digits = 3), "Calibrated (100 sim)"),
    paste(round(res_ss_ew_150sim$AUC, digits = 3), "Calibrated (150 sim)"),
    paste(round(res_ss_ew_200sim$AUC, digits = 3), "Calibrated (200 sim)"),
    paste(round(res_random_ew$AUC, digits = 3), "Random")), cex = 0.9)
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(x = res_ss_mw_50sim$roc_stats$FPR, y = res_ss_mw_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Model-wise synergies (HSA)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_ss_mw_100sim$roc_stats$FPR, y = res_ss_mw_100sim$roc_stats$TPR, 
  lwd = 3, col = my_palette[2])
lines(x = res_ss_mw_150sim$roc_stats$FPR, y = res_ss_mw_150sim$roc_stats$TPR, 
  lwd = 3, col = my_palette[3])
lines(x = res_ss_mw_200sim$roc_stats$FPR, y = res_ss_mw_200sim$roc_stats$TPR, 
  lwd = 3, col = my_palette[4])
lines(x = res_random_mw$roc_stats$FPR, y = res_random_mw$roc_stats$TPR, 
  lwd = 3, col = my_palette[5])
legend('bottomright', title = 'AUC', col = my_palette[1:5], pch = 19,
  legend = c(paste(round(res_ss_mw_50sim$AUC, digits = 3), "Calibrated (50 sim)"),
    paste(round(res_ss_mw_100sim$AUC, digits = 3), "Calibrated (100 sim)"),
    paste(round(res_ss_mw_150sim$AUC, digits = 3), "Calibrated (150 sim)"),
    paste(round(res_ss_mw_200sim$AUC, digits = 3), "Calibrated (200 sim)"),
    paste(round(res_random_mw$AUC, digits = 3), "Random")), cex = 0.9)
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
```

<img src="index_files/figure-html/ROC HSA Cascade 2.0-1.png" width="50%" /><img src="index_files/figure-html/ROC HSA Cascade 2.0-2.png" width="50%" />

### PR curves {-}


```r
pr_ss_ew_hsa_50sim = pr.curve(scores.class0 = pred_ew_hsa %>% pull(ss_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_ss_ew_hsa_100sim = pr.curve(scores.class0 = pred_ew_hsa %>% pull(ss_score_100sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_hsa %>% pull(observed), curve = TRUE)
pr_ss_ew_hsa_150sim = pr.curve(scores.class0 = pred_ew_hsa %>% pull(ss_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_hsa %>% pull(observed), curve = TRUE)
pr_ss_ew_hsa_200sim = pr.curve(scores.class0 = pred_ew_hsa %>% pull(ss_score_200sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_hsa %>% pull(observed), curve = TRUE)
pr_random_ew_hsa = pr.curve(scores.class0 = pred_ew_hsa %>% pull(random_score) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_hsa %>% pull(observed), curve = TRUE)

pr_ss_mw_hsa_50sim = pr.curve(scores.class0 = pred_mw_hsa %>% pull(synergy_prob_ss_50sim), 
  weights.class0 = pred_mw_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_ss_mw_hsa_100sim = pr.curve(scores.class0 = pred_mw_hsa %>% pull(synergy_prob_ss_100sim), 
  weights.class0 = pred_mw_hsa %>% pull(observed), curve = TRUE)
pr_ss_mw_hsa_150sim = pr.curve(scores.class0 = pred_mw_hsa %>% pull(synergy_prob_ss_150sim), 
  weights.class0 = pred_mw_hsa %>% pull(observed), curve = TRUE)
pr_ss_mw_hsa_200sim = pr.curve(scores.class0 = pred_mw_hsa %>% pull(synergy_prob_ss_200sim), 
  weights.class0 = pred_mw_hsa %>% pull(observed), curve = TRUE)
pr_random_mw_hsa = pr.curve(scores.class0 = pred_mw_hsa %>% pull(synergy_prob_random), 
  weights.class0 = pred_mw_hsa %>% pull(observed), curve = TRUE)

plot(pr_ss_ew_hsa_50sim, main = 'PR curve, Ensemble-wise synergies (HSA)', 
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_ss_ew_hsa_100sim, add = TRUE, color = my_palette[2])
plot(pr_ss_ew_hsa_150sim, add = TRUE, color = my_palette[3])
plot(pr_ss_ew_hsa_200sim, add = TRUE, color = my_palette[4])
plot(pr_random_ew_hsa, add = TRUE, color = my_palette[5])
legend('topright', title = 'AUC', col = my_palette[1:5], pch = 19,
  legend = c(paste(round(pr_ss_ew_hsa_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_ss_ew_hsa_100sim$auc.davis.goadrich, digits = 3), "Calibrated (100 sim)"),
    paste(round(pr_ss_ew_hsa_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_ss_ew_hsa_200sim$auc.davis.goadrich, digits = 3), "Calibrated (200 sim)"),
    paste(round(pr_random_ew_hsa$auc.davis.goadrich, digits = 3), "Random")))
grid(lwd = 0.5)

plot(pr_ss_mw_hsa_50sim, main = 'PR curve, Model-wise synergies (HSA)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_ss_mw_hsa_100sim, add = TRUE, color = my_palette[2])
plot(pr_ss_mw_hsa_150sim, add = TRUE, color = my_palette[3])
plot(pr_ss_mw_hsa_200sim, add = TRUE, color = my_palette[4])
plot(pr_random_mw_hsa, add = TRUE, color = my_palette[5])
legend('topright', title = 'AUC', col = my_palette[1:5], pch = 19,
  legend = c(paste(round(pr_ss_mw_hsa_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_ss_mw_hsa_100sim$auc.davis.goadrich, digits = 3), "Calibrated (100 sim)"),
    paste(round(pr_ss_mw_hsa_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_ss_mw_hsa_200sim$auc.davis.goadrich, digits = 3), "Calibrated (200 sim)"),
    paste(round(pr_random_mw_hsa$auc.davis.goadrich, digits = 3), "Random")))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/PR HSA Cascade 2.0-1.png" width="50%" /><img src="index_files/figure-html/PR HSA Cascade 2.0-2.png" width="50%" />

When I saw the above I was like:  
<iframe src="https://giphy.com/embed/l3q2K5jinAlChoCLS" width="250" height="270" frameBorder="0" class="giphy-embed" allowFullScreen></iframe>

:::{.green-box}
- PR curves show the **true colors** of our models' prediction performance for our imbalanced dataset!
- Both ROC & PR curves show that **random models have pretty much equal performance to calibrated models**
- The performance *scaling* is still consistent with the PR AUCs but on a much smaller scale! 
- The **model-wise approach is slightly better than the ensemble-wise**
:::

### ROC AUC sensitivity {-}

Investigate same thing as described in [here](#auc-sensitivity).
We will combine the synergy scores from the random simulations with the results from the $150$ Gitsbe simulations.


```r
# Ensemble-wise
betas = seq(from = -10, to = 10, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_ew_hsa = pred_ew_hsa %>% mutate(combined_score = ss_score_150sim + beta * random_score)
  res = get_roc_stats(df = pred_ew_hsa, pred_col = "combined_score", label_col = "observed")
  auc_value = res$AUC
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2]) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/AUC sensitivity (HSA - cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

auc_values_mw = sapply(weights, function(w) {
  pred_mw_hsa = pred_mw_hsa %>% 
    mutate(weighted_prob = (1 - w) * pred_mw_hsa$synergy_prob_ss_150sim + w * pred_mw_hsa$synergy_prob_random)
  res = get_roc_stats(df = pred_mw_hsa, pred_col = "weighted_prob", label_col = "observed", direction = ">")
  auc_value = res$AUC
})

df_mw = as_tibble(cbind(weights, auc_values_mw))

ggline(data = df_mw, x = "weights", y = "auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids() #+ ylim(0.4, 0.9)
```

<img src="index_files/figure-html/AUC sensitivity (HSA - cascade 2.0)-2.png" width="80%" style="display: block; margin: auto;" />

:::{.green-box}
- Symmetricity (Ensemble-wise): $AUC_{\beta \rightarrow +\infty} + AUC_{\beta \rightarrow -\infty} \approx 1$
- Random models perform worse than calibrated ones (though difference is very small)
- There are $\beta$ values that can boost the predictive performance of the combined synergy classifier and a $w$ weight in the model-wise case (though the significance in performance gain is negligible).
:::

### PR AUC sensitivity {-}

Again we try a linear combination of the predicted scores from the random simulations with the results from the $150$ Gitsbe simulations.


```r
# Ensemble-wise
betas = seq(from = -20, to = 20, by = 0.1)

pr_auc_values_ew = sapply(betas, function(beta) {
  pred_ew_hsa = pred_ew_hsa %>% mutate(combined_score = ss_score_150sim + beta * random_score)
  res = pr.curve(scores.class0 = pred_ew_hsa %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_ew_hsa %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_pr_ew = as_tibble(cbind(betas, pr_auc_values_ew))

ggline(data = df_pr_ew, x = "betas", y = "pr_auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "PR-AUC (Area Under PR Curve)",
  title = TeX("PR-AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2]) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/PR AUC sensitivity (HSA - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

pr_auc_values_mw = sapply(weights, function(w) {
  pred_mw_hsa = pred_mw_hsa %>% 
    mutate(weighted_prob = (1 - w) * pred_mw_hsa$synergy_prob_ss_150sim + w * pred_mw_hsa$synergy_prob_random)
  res = pr.curve(scores.class0 = pred_mw_hsa %>% pull(weighted_prob), 
    weights.class0 = pred_mw_hsa %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_pr_mw = as_tibble(cbind(weights, pr_auc_values_mw))

ggline(data = df_pr_mw, x = "weights", y = "pr_auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "PR-AUC (Area Under PR Curve)",
  title = TeX("PR-AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/PR AUC sensitivity (HSA - Cascade 2.0)-2.png" width="80%" style="display: block; margin: auto;" />

### Logistic Regression {-}

We tried fitting a model using logistic regression as a different approach to augment the results from calibrated simulations with the random ones (**HSA, ensemble-wise results**).


```r
model = glm(formula = observed ~ ss_score_150sim + random_score - 1, data = pred_ew_hsa, family = binomial())
model_tidy = broom::tidy(model)
coef1 = model_tidy %>% filter(term == "ss_score_150sim") %>% pull(estimate)
coef2 = model_tidy %>% filter(term == "random_score") %>% pull(estimate)
pred_ew_hsa = pred_ew_hsa %>% mutate(glm = coef1 * ss_score_150sim + coef2 * random_score)
res_roc = get_roc_stats(df = pred_ew_hsa, pred_col = "glm", label_col = "observed")
res_pr = pr.curve(scores.class0 = pred_ew_hsa %>% pull(glm) %>% (function(x) {-x}), weights.class0 = pred_ew_hsa %>% pull(observed))
```

The model with the coefficients is as follows (note that adding an intercept makes ROC AUC result considerably worse):

```r
extract_eq(model, use_coefs = TRUE)
```

$$
\log\left[ \frac { P( \text{observed} = \text{1} ) }{ 1 - P( \text{observed} = \text{1} ) } \right] = 13.37(\text{ss\_score\_150sim}) + 11.83(\text{random\_score}) + \epsilon
$$

:::{.orange-box}
The ROC AUC produced with a logistic regression model is lower than the calibrated models (with $150$ Gitsbe simulations): 0.6938776 (PR-AUC is also lower: 0.0652647).
:::

### Regularized Logistic Regression {-}

We try a regularized logistic regression approach using the `glmnet` R package [@R-glmnet].
We cross validate the $\lambda$ parameter and try with different $\alpha \in [0,1]$ ($\alpha=0$ means Ridge regression, $\alpha=1$ means LASSO, in between means Elastic net) while either minimizing the missclassification error (`type.measure="class"`) or maximizing the ROC-AUC (`type.measure = "auc"`).
For each respective $\alpha$ we choose the $\lambda_{min}$ as the one the minimizes the average CV error.
The intercept was again excluded as it resulted in worse AUC performance.


```r
x = pred_ew_hsa %>% select(ss_score_150sim, random_score) %>% as.matrix()
y = pred_ew_hsa %>% pull(observed)

data_list = list()
index = 1
for (i in 0:10) { # from Ridge to LASSO
  a = i/10
  for (measure in c("auc", "class")) {
    set.seed(42) # for reproducibility
    cvfit = cv.glmnet(x, y, family = "binomial", type.measure = measure, intercept = FALSE, alpha = a)
    coef_mat = coef(cvfit, s = "lambda.min")
    pred_ew_hsa = pred_ew_hsa %>% mutate(glm_reg = coef_mat[1] + coef_mat[2] * ss_score_150sim + coef_mat[3] * random_score)
    res_roc = get_roc_stats(df = pred_ew_hsa, pred_col = "glm_reg", label_col = "observed")
    pr_roc = pr.curve(scores.class0 = pred_ew_hsa %>% pull(glm_reg) %>% (function(x) {-x}), 
      weights.class0 = pred_ew_hsa %>% pull(observed))
    data_list[[index]] = as_tibble_row(list(alpha = a, measure = measure, ROC_AUC = res_roc$AUC, PR_AUC = pr_roc$auc.davis.goadrich))
    index = index + 1
  }
}

data = bind_rows(data_list)

# List the best two results
data %>% arrange(desc(ROC_AUC)) %>% slice(1:4) %>% kable()
```



 alpha  measure      ROC_AUC      PR_AUC
------  --------  ----------  ----------
   0.0  class      0.6870748   0.0634616
   0.0  auc        0.6848073   0.0629667
   0.1  auc        0.6836735   0.0625847
   0.5  auc        0.6768707   0.0655809

:::{.orange-box}
The best ROC AUC produced with a regularized logistic regression model is also lower than the one using calibrated models alone (with $150$ Gitsbe simulations).

Note that we get warnings when using `glmnet` because of the small number of observations for the positive class (observed synergies).
Resulting coefficients vary, but tend to be either all too small or **larger on the random model predictor**.
:::

### MAMSE ROC {-}

Using the `MAMSE` R package [@R-MAMSE] we try another method to combine the predictor values from the calibrated and random models.
The resulting ROC curve gets a little bit distored and AUC is not statistically better from the reference sample population (i.e. the calibrated `Gitsbe` models with $150$):


```r
# healthy => non-synergy, diseased => synergy
healthy = list()
healthy[[1]] = pred_ew_hsa %>% filter(observed == 0) %>% pull(ss_score_150sim)
healthy[[2]] = pred_ew_hsa %>% filter(observed == 0) %>% pull(random_score)

diseased = list()
diseased[[1]] = pred_ew_hsa %>% filter(observed == 1) %>% pull(ss_score_150sim)
diseased[[2]] = pred_ew_hsa %>% filter(observed == 1) %>% pull(random_score)

plot(roc(healthy = healthy, diseased = diseased, smalldiseased=TRUE, AUC=TRUE,
  wh=NULL, wd=NULL, FPR=NULL, method="np"))
```

<img src="index_files/figure-html/MAMSE test-1.png" width="672" />

## Calibrated vs Random (Bliss) {-}

:::{.note}
- *Bliss* refers to the synergy method used in `Drabme` to assess the synergies from the `gitsbe` models
- We test performance using ROC AUC for both the *ensemble-wise* and *model-wise* synergies from `Drabme`
- **Calibrated** models: fitted to steady state
- **Random** models: produced via `abmlog` (see [here](#random-model-results) and used in `Drabme` with `synergy_method: bliss`
:::

Load results:

```r
# 'ss' => calibrated models, 'random' => random models

## Bliss results
ss_bliss_ensemblewise_10sim_file = paste0("results/bliss/cascade_2.0_ss_10sim_fixpoints_ensemblewise_synergies.tab")
ss_bliss_modelwise_10sim_file = paste0("results/bliss/cascade_2.0_ss_10sim_fixpoints_modelwise_synergies.tab")
ss_bliss_ensemblewise_30sim_file = paste0("results/bliss/cascade_2.0_ss_30sim_fixpoints_ensemblewise_synergies.tab")
ss_bliss_modelwise_30sim_file = paste0("results/bliss/cascade_2.0_ss_30sim_fixpoints_modelwise_synergies.tab")
ss_bliss_ensemblewise_50sim_file = paste0("results/bliss/cascade_2.0_ss_50sim_fixpoints_ensemblewise_synergies.tab")
ss_bliss_modelwise_50sim_file = paste0("results/bliss/cascade_2.0_ss_50sim_fixpoints_modelwise_synergies.tab")
ss_bliss_ensemblewise_70sim_file = paste0("results/bliss/cascade_2.0_ss_70sim_fixpoints_ensemblewise_synergies.tab")
ss_bliss_modelwise_70sim_file = paste0("results/bliss/cascade_2.0_ss_70sim_fixpoints_modelwise_synergies.tab")
ss_bliss_ensemblewise_100sim_file = paste0("results/bliss/cascade_2.0_ss_100sim_fixpoints_ensemblewise_synergies.tab")
ss_bliss_modelwise_100sim_file = paste0("results/bliss/cascade_2.0_ss_100sim_fixpoints_modelwise_synergies.tab")
ss_bliss_ensemblewise_150sim_file = paste0("results/bliss/cascade_2.0_ss_150sim_fixpoints_ensemblewise_synergies.tab")
ss_bliss_modelwise_150sim_file = paste0("results/bliss/cascade_2.0_ss_150sim_fixpoints_modelwise_synergies.tab")
ss_bliss_ensemblewise_200sim_file = paste0("results/bliss/cascade_2.0_ss_200sim_fixpoints_ensemblewise_synergies.tab")
ss_bliss_modelwise_200sim_file = paste0("results/bliss/cascade_2.0_ss_200sim_fixpoints_modelwise_synergies.tab")

random_bliss_ensemblewise_file = paste0("results/bliss/cascade_2.0_random_bliss_ensemblewise_synergies.tab")
random_bliss_modelwise_file = paste0("results/bliss/cascade_2.0_random_bliss_modelwise_synergies.tab")

ss_bliss_ensemblewise_synergies_10sim = emba::get_synergy_scores(ss_bliss_ensemblewise_10sim_file)
ss_bliss_modelwise_synergies_10sim = emba::get_synergy_scores(ss_bliss_modelwise_10sim_file, file_type = "modelwise")
ss_bliss_ensemblewise_synergies_30sim = emba::get_synergy_scores(ss_bliss_ensemblewise_30sim_file)
ss_bliss_modelwise_synergies_30sim = emba::get_synergy_scores(ss_bliss_modelwise_30sim_file, file_type = "modelwise")
ss_bliss_ensemblewise_synergies_50sim = emba::get_synergy_scores(ss_bliss_ensemblewise_50sim_file)
ss_bliss_modelwise_synergies_50sim = emba::get_synergy_scores(ss_bliss_modelwise_50sim_file, file_type = "modelwise")
ss_bliss_ensemblewise_synergies_70sim = emba::get_synergy_scores(ss_bliss_ensemblewise_70sim_file)
ss_bliss_modelwise_synergies_70sim = emba::get_synergy_scores(ss_bliss_modelwise_70sim_file, file_type = "modelwise")
ss_bliss_ensemblewise_synergies_100sim = emba::get_synergy_scores(ss_bliss_ensemblewise_100sim_file)
ss_bliss_modelwise_synergies_100sim = emba::get_synergy_scores(ss_bliss_modelwise_100sim_file, file_type = "modelwise")
ss_bliss_ensemblewise_synergies_150sim = emba::get_synergy_scores(ss_bliss_ensemblewise_150sim_file)
ss_bliss_modelwise_synergies_150sim = emba::get_synergy_scores(ss_bliss_modelwise_150sim_file, file_type = "modelwise")
ss_bliss_ensemblewise_synergies_200sim = emba::get_synergy_scores(ss_bliss_ensemblewise_200sim_file)
ss_bliss_modelwise_synergies_200sim = emba::get_synergy_scores(ss_bliss_modelwise_200sim_file, file_type = "modelwise")

random_bliss_ensemblewise_synergies = emba::get_synergy_scores(random_bliss_ensemblewise_file)
random_bliss_modelwise_synergies = emba::get_synergy_scores(random_bliss_modelwise_file, file_type = "modelwise")

# calculate probability of synergy in the modelwise results
ss_bliss_modelwise_synergies_10sim = ss_bliss_modelwise_synergies_10sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
ss_bliss_modelwise_synergies_30sim = ss_bliss_modelwise_synergies_30sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
ss_bliss_modelwise_synergies_50sim = ss_bliss_modelwise_synergies_50sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
ss_bliss_modelwise_synergies_70sim = ss_bliss_modelwise_synergies_70sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
ss_bliss_modelwise_synergies_100sim = ss_bliss_modelwise_synergies_100sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
ss_bliss_modelwise_synergies_150sim = ss_bliss_modelwise_synergies_150sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
ss_bliss_modelwise_synergies_200sim = ss_bliss_modelwise_synergies_200sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
random_bliss_modelwise_synergies = random_bliss_modelwise_synergies %>%
  mutate(synergy_prob_random = synergies/(synergies + `non-synergies`))
```

### ROC curves {-}


```r
# 'ew' => ensemble-wise, 'mw' => model-wise
pred_ew_bliss = bind_cols(ss_bliss_ensemblewise_synergies_10sim %>% rename(ss_score_10sim = score), 
  ss_bliss_ensemblewise_synergies_30sim %>% select(score) %>% rename(ss_score_30sim = score),
  ss_bliss_ensemblewise_synergies_50sim %>% select(score) %>% rename(ss_score_50sim = score),
  ss_bliss_ensemblewise_synergies_70sim %>% select(score) %>% rename(ss_score_70sim = score),
  ss_bliss_ensemblewise_synergies_100sim %>% select(score) %>% rename(ss_score_100sim = score),
  ss_bliss_ensemblewise_synergies_150sim %>% select(score) %>% rename(ss_score_150sim = score),
  ss_bliss_ensemblewise_synergies_200sim %>% select(score) %>% rename(ss_score_200sim = score),
  random_bliss_ensemblewise_synergies %>% select(score) %>% rename(random_score = score), 
  as_tibble_col(observed, column_name = "observed"))

pred_mw_bliss = bind_cols(
  ss_bliss_modelwise_synergies_10sim %>% select(perturbation, synergy_prob_ss) %>% rename(synergy_prob_ss_10sim = synergy_prob_ss),
  ss_bliss_modelwise_synergies_30sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_30sim = synergy_prob_ss),
  ss_bliss_modelwise_synergies_50sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_50sim = synergy_prob_ss),
  ss_bliss_modelwise_synergies_70sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_70sim = synergy_prob_ss),
  ss_bliss_modelwise_synergies_100sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_100sim = synergy_prob_ss),
  ss_bliss_modelwise_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_150sim = synergy_prob_ss),
  ss_bliss_modelwise_synergies_200sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_200sim = synergy_prob_ss),
  random_bliss_modelwise_synergies %>% select(synergy_prob_random),
  as_tibble_col(observed, column_name = "observed"))

res_ss_ew_10sim = get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score_10sim", label_col = "observed")
res_ss_ew_30sim = get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score_30sim", label_col = "observed")
res_ss_ew_50sim = get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score_50sim", label_col = "observed")
res_ss_ew_70sim = get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score_70sim", label_col = "observed")
res_ss_ew_100sim = get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score_100sim", label_col = "observed")
res_ss_ew_150sim = get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score_150sim", label_col = "observed")
res_ss_ew_200sim = get_roc_stats(df = pred_ew_bliss, pred_col = "ss_score_200sim", label_col = "observed")
res_random_ew = get_roc_stats(df = pred_ew_bliss, pred_col = "random_score", label_col = "observed")

res_ss_mw_10sim = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_ss_10sim", label_col = "observed", direction = ">")
res_ss_mw_30sim = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_ss_30sim", label_col = "observed", direction = ">")
res_ss_mw_50sim = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_ss_50sim", label_col = "observed", direction = ">")
res_ss_mw_70sim = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_ss_70sim", label_col = "observed", direction = ">")
res_ss_mw_100sim = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_ss_100sim", label_col = "observed", direction = ">")
res_ss_mw_150sim = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_ss_150sim", label_col = "observed", direction = ">")
res_ss_mw_200sim = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_ss_200sim", label_col = "observed", direction = ">")
res_random_mw = get_roc_stats(df = pred_mw_bliss, pred_col = "synergy_prob_random", label_col = "observed", direction = ">")

# Plot ROCs
plot(x = res_ss_ew_10sim$roc_stats$FPR, y = res_ss_ew_10sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_ss_ew_30sim$roc_stats$FPR, y = res_ss_ew_30sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = res_ss_ew_50sim$roc_stats$FPR, y = res_ss_ew_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = res_ss_ew_70sim$roc_stats$FPR, y = res_ss_ew_70sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
lines(x = res_ss_ew_100sim$roc_stats$FPR, y = res_ss_ew_100sim$roc_stats$TPR,
  lwd = 3, col = my_palette[5])
lines(x = res_ss_ew_150sim$roc_stats$FPR, y = res_ss_ew_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[6])
lines(x = res_ss_ew_200sim$roc_stats$FPR, y = res_ss_ew_200sim$roc_stats$TPR,
  lwd = 3, col = my_palette[7])
lines(x = res_random_ew$roc_stats$FPR, y = res_random_ew$roc_stats$TPR, 
  lwd = 3, col = my_palette[8])
legend('bottomright', title = 'AUC', col = my_palette[1:8], pch = 19,
  legend = c(paste(round(res_ss_ew_10sim$AUC, digits = 2), "Calibrated (10 sim)"),
    paste(round(res_ss_ew_30sim$AUC, digits = 2), "Calibrated (30 sim)"),
    paste(round(res_ss_ew_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(res_ss_ew_70sim$AUC, digits = 2), "Calibrated (70 sim)"),
    paste(round(res_ss_ew_100sim$AUC, digits = 2), "Calibrated (100 sim)"),
    paste(round(res_ss_ew_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(res_ss_ew_200sim$AUC, digits = 2), "Calibrated (200 sim)"),
    paste(round(res_random_ew$AUC, digits = 2), "Random")), cex = 0.7)
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(x = res_ss_mw_10sim$roc_stats$FPR, y = res_ss_mw_10sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Model-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = res_ss_mw_30sim$roc_stats$FPR, y = res_ss_mw_30sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = res_ss_mw_50sim$roc_stats$FPR, y = res_ss_mw_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = res_ss_mw_70sim$roc_stats$FPR, y = res_ss_mw_70sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
lines(x = res_ss_mw_100sim$roc_stats$FPR, y = res_ss_mw_100sim$roc_stats$TPR,
  lwd = 3, col = my_palette[5])
lines(x = res_ss_mw_150sim$roc_stats$FPR, y = res_ss_mw_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[6])
lines(x = res_ss_mw_200sim$roc_stats$FPR, y = res_ss_mw_200sim$roc_stats$TPR,
  lwd = 3, col = my_palette[7])
lines(x = res_random_mw$roc_stats$FPR, y = res_random_mw$roc_stats$TPR, 
  lwd = 3, col = my_palette[8])
legend('bottomright', title = 'AUC', col = my_palette[1:8], pch = 19,
  legend = c(paste(round(res_ss_mw_10sim$AUC, digits = 2), "Calibrated (10 sim)"),
    paste(round(res_ss_mw_30sim$AUC, digits = 2), "Calibrated (30 sim)"),
    paste(round(res_ss_mw_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(res_ss_mw_70sim$AUC, digits = 2), "Calibrated (70 sim)"),
    paste(round(res_ss_mw_100sim$AUC, digits = 2), "Calibrated (100 sim)"),
    paste(round(res_ss_mw_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(res_ss_mw_200sim$AUC, digits = 2), "Calibrated (200 sim)"),
    paste(round(res_random_mw$AUC, digits = 2), "Random")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
```

<img src="index_files/figure-html/ROC Bliss Cascade 2.0-1.png" width="50%" /><img src="index_files/figure-html/ROC Bliss Cascade 2.0-2.png" width="50%" />

:::{.green-box}
- Model-wise results *scale* with respect to the number of `Gitsbe` simulations (more **calibrated** models, better performance).
- Ensemble-wise performance is disproportionate compared to the model-wise performance 
:::

### PR curves {-}


```r
pr_ss_ew_bliss_10sim = pr.curve(scores.class0 = pred_ew_bliss %>% pull(ss_score_10sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_ss_ew_bliss_30sim = pr.curve(scores.class0 = pred_ew_bliss %>% pull(ss_score_30sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)
pr_ss_ew_bliss_50sim = pr.curve(scores.class0 = pred_ew_bliss %>% pull(ss_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)
pr_ss_ew_bliss_70sim = pr.curve(scores.class0 = pred_ew_bliss %>% pull(ss_score_70sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)
pr_ss_ew_bliss_100sim = pr.curve(scores.class0 = pred_ew_bliss %>% pull(ss_score_100sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)
pr_ss_ew_bliss_150sim = pr.curve(scores.class0 = pred_ew_bliss %>% pull(ss_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)
pr_ss_ew_bliss_200sim = pr.curve(scores.class0 = pred_ew_bliss %>% pull(ss_score_200sim) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)
pr_random_ew_bliss = pr.curve(scores.class0 = pred_ew_bliss %>% pull(random_score) %>% (function(x) {-x}), 
  weights.class0 = pred_ew_bliss %>% pull(observed), curve = TRUE)

pr_ss_mw_bliss_10sim = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_ss_10sim), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_ss_mw_bliss_30sim = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_ss_30sim), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE)
pr_ss_mw_bliss_50sim = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_ss_50sim), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE)
pr_ss_mw_bliss_70sim = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_ss_70sim), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE)
pr_ss_mw_bliss_100sim = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_ss_100sim), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE)
pr_ss_mw_bliss_150sim = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_ss_150sim), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE)
pr_ss_mw_bliss_200sim = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_ss_200sim), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE)
pr_random_mw_bliss = pr.curve(scores.class0 = pred_mw_bliss %>% pull(synergy_prob_random), 
  weights.class0 = pred_mw_bliss %>% pull(observed), curve = TRUE)

plot(pr_ss_ew_bliss_10sim, main = 'PR curve, Ensemble-wise synergies (Bliss)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_ss_ew_bliss_30sim, add = TRUE, color = my_palette[2])
plot(pr_ss_ew_bliss_50sim, add = TRUE, color = my_palette[3])
plot(pr_ss_ew_bliss_70sim, add = TRUE, color = my_palette[4])
plot(pr_ss_ew_bliss_100sim, add = TRUE, color = my_palette[5])
plot(pr_ss_ew_bliss_150sim, add = TRUE, color = my_palette[6])
plot(pr_ss_ew_bliss_200sim, add = TRUE, color = my_palette[7])
plot(pr_random_ew_bliss, add = TRUE, color = my_palette[8])
legend('topright', title = 'AUC', col = my_palette[1:8], pch = 19,
  legend = c(paste(round(pr_ss_ew_bliss_10sim$auc.davis.goadrich, digits = 3), "Calibrated (10 sim)"),
    paste(round(pr_ss_ew_bliss_30sim$auc.davis.goadrich, digits = 3), "Calibrated (30 sim)"),
    paste(round(pr_ss_ew_bliss_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_ss_ew_bliss_70sim$auc.davis.goadrich, digits = 3), "Calibrated (70 sim)"),
    paste(round(pr_ss_ew_bliss_100sim$auc.davis.goadrich, digits = 3), "Calibrated (100 sim)"),
    paste(round(pr_ss_ew_bliss_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_ss_ew_bliss_200sim$auc.davis.goadrich, digits = 3), "Calibrated (200 sim)"),
    paste(round(pr_random_ew_bliss$auc.davis.goadrich, digits = 3), "Random")))
grid(lwd = 0.5)

plot(pr_ss_mw_bliss_10sim, main = 'PR curve, Model-wise synergies (Bliss)', 
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_ss_mw_bliss_30sim, add = TRUE, color = my_palette[2])
plot(pr_ss_mw_bliss_50sim, add = TRUE, color = my_palette[3])
plot(pr_ss_mw_bliss_70sim, add = TRUE, color = my_palette[4])
plot(pr_ss_mw_bliss_100sim, add = TRUE, color = my_palette[5])
plot(pr_ss_mw_bliss_150sim, add = TRUE, color = my_palette[6])
plot(pr_ss_mw_bliss_200sim, add = TRUE, color = my_palette[7])
plot(pr_random_mw_bliss, add = TRUE, color = my_palette[8])
legend('topright', title = 'AUC', col = my_palette[1:8], pch = 19,
  legend = c(paste(round(pr_ss_mw_bliss_10sim$auc.davis.goadrich, digits = 3), "Calibrated (10 sim)"),
    paste(round(pr_ss_mw_bliss_30sim$auc.davis.goadrich, digits = 3), "Calibrated (30 sim)"),
    paste(round(pr_ss_mw_bliss_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_ss_mw_bliss_70sim$auc.davis.goadrich, digits = 3), "Calibrated (70 sim)"),
    paste(round(pr_ss_mw_bliss_100sim$auc.davis.goadrich, digits = 3), "Calibrated (100 sim)"),
    paste(round(pr_ss_mw_bliss_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_ss_mw_bliss_200sim$auc.davis.goadrich, digits = 3), "Calibrated (200 sim)"),
    paste(round(pr_random_mw_bliss$auc.davis.goadrich, digits = 3), "Random")))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/PR Bliss Cascade 2.0-1.png" width="50%" /><img src="index_files/figure-html/PR Bliss Cascade 2.0-2.png" width="50%" />

### ROC AUC sensitivity {-}

Investigate same thing as described in [here](#auc-sensitivity).
We will combine the synergy scores from the random simulations with the results from the $50$ Gitsbe simulations.


```r
# Ensemble-wise
betas = seq(from = -20, to = 20, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_ew_bliss = pred_ew_bliss %>% mutate(combined_score = ss_score_50sim + beta * random_score)
  res = get_roc_stats(df = pred_ew_bliss, pred_col = "combined_score", label_col = "observed")
  auc_value = res$AUC
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2]) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/AUC sensitivity (Bliss - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

auc_values_mw = sapply(weights, function(w) {
  pred_mw_bliss = pred_mw_bliss %>% 
    mutate(weighted_prob = (1 - w) * pred_mw_bliss$synergy_prob_ss_50sim + w * pred_mw_bliss$synergy_prob_random)
  res = get_roc_stats(df = pred_mw_bliss, pred_col = "weighted_prob", label_col = "observed", direction = ">")
  auc_value = res$AUC
})

df_mw = as_tibble(cbind(weights, auc_values_mw))

ggline(data = df_mw, x = "weights", y = "auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/AUC sensitivity (Bliss - Cascade 2.0)-2.png" width="80%" style="display: block; margin: auto;" />

:::{.green-box}
- Symmetricity (Ensemble-wise): $AUC_{\beta \rightarrow +\infty} + AUC_{\beta \rightarrow -\infty} \approx 1$
- **Random models perform better than calibrated ones**
- Combining the synergy results using the weighted probability score does not bring any significant difference in performance
- Using the $\beta$ parameter to boost the ensemble synergy results does not work for the Bliss results (improvement was seen only for the HSA results)
:::

### PR AUC sensitivity {-}

Again we try a linear combination of the predicted scores from the random simulations with the results from the $150$ Gitsbe simulations.


```r
# Ensemble-wise
betas = seq(from = -20, to = 20, by = 0.1)

pr_auc_values_ew = sapply(betas, function(beta) {
  pred_ew_bliss = pred_ew_bliss %>% mutate(combined_score = ss_score_150sim + beta * random_score)
  res = pr.curve(scores.class0 = pred_ew_bliss %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_ew_bliss %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_pr_ew = as_tibble(cbind(betas, pr_auc_values_ew))

ggline(data = df_pr_ew, x = "betas", y = "pr_auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "PR-AUC (Area Under PR Curve)",
  title = TeX("PR-AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2]) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/PR AUC sensitivity (Bliss - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

pr_auc_values_mw = sapply(weights, function(w) {
  pred_mw_bliss = pred_mw_bliss %>% 
    mutate(weighted_prob = (1 - w) * pred_mw_bliss$synergy_prob_ss_150sim + w * pred_mw_bliss$synergy_prob_random)
  res = pr.curve(scores.class0 = pred_mw_bliss %>% pull(weighted_prob), 
    weights.class0 = pred_mw_bliss %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_pr_mw = as_tibble(cbind(weights, pr_auc_values_mw))

ggline(data = df_pr_mw, x = "weights", y = "pr_auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "PR-AUC (Area Under PR Curve)",
  title = TeX("PR-AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/PR AUC sensitivity (Bliss - Cascade 2.0)-2.png" width="80%" style="display: block; margin: auto;" />

## Correlation {-}

We test for correlation between some of the results shown in the ROC curves (the calibrated models are from the $150$ simulation results).
This means *ensemble-wise* vs *model-wise*, *random* models vs *calibrated (ss)* models and *HSA* vs *Bliss* synergy assessment.
*P-values* are represented at 3 significant levels: $0.05, 0.01, 0.001$ (\*, \*\*, \*\*\*) and the correlation coefficient is calculated using Kendall's *tau* statistic.


```r
synergy_scores = bind_cols(
  pred_ew_hsa %>% select(ss_score_100sim, random_score) %>% rename(ss_ensemble_hsa = ss_score_100sim, random_ensemble_hsa = random_score),
  pred_ew_bliss %>% select(ss_score_100sim, random_score) %>% rename(ss_ensemble_bliss = ss_score_100sim, random_ensemble_bliss = random_score),
  pred_mw_hsa %>% select(synergy_prob_ss_100sim, synergy_prob_random) %>% 
    rename(ss_modelwise_hsa = synergy_prob_ss_100sim, random_modelwise_hsa = synergy_prob_random),
  pred_mw_bliss %>% select(synergy_prob_ss_100sim, synergy_prob_random) %>% 
    rename(ss_modelwise_bliss = synergy_prob_ss_100sim, random_modelwise_bliss = synergy_prob_random)
  )

M = cor(synergy_scores, method = "kendall")
res = cor.mtest(mat = synergy_scores, method = "kendall")
corrplot(corr = M, type ="upper", p.mat = res$p, sig.level = c(.001, .01, .05), 
  pch.cex = 1, pch.col = "white", insig = "label_sig", tl.col = "black", tl.srt = 45)
```

<img src="index_files/figure-html/Correlation of ROC results (Cascade 2.0)-1.png" width="2100" />

:::{.green-box}
- **HSA and Bliss results correlate**, especially for the model-wise results (ensemble-wise correlation is not as strong).
- **Model-wise don't correlate with ensemble-wise results**.
:::

## Topology Mutations {-}

:::{.note}
We run `Gitsbe` simulations with $50$ topology mutations (bootstrap value, reduced to $10$ after models with stabla states have been found), both for $50$ and $150$ simulations and both **fitting to steady state** (calibrated models) and to a **proliferative phenotype** (so not random models but as close as it can get to that since we are discussing topology mutations).

Results were assessed via `Drabme`'s ensemble and modelwise results, using both the HSA and Bliss method.
:::

Load the results:


```r
# 'ss' => calibrated models, 'rand' => proliferative models (so not random but kind of!)
# 'ew' => ensemble-wise, 'mw' => modelwise

## HSA results ss
topo_ss_hsa_ew_50sim_file = paste0("results/topology-only/cascade_2.0_ss_50sim_fixpoints_hsa_ensemblewise_synergies.tab")
topo_ss_hsa_mw_50sim_file = paste0("results/topology-only/cascade_2.0_ss_50sim_fixpoints_hsa_modelwise_synergies.tab")
topo_ss_hsa_ew_150sim_file = paste0("results/topology-only/cascade_2.0_ss_150sim_fixpoints_hsa_ensemblewise_synergies.tab")
topo_ss_hsa_mw_150sim_file = paste0("results/topology-only/cascade_2.0_ss_150sim_fixpoints_hsa_modelwise_synergies.tab")

topo_ss_hsa_ew_synergies_50sim = emba::get_synergy_scores(topo_ss_hsa_ew_50sim_file)
topo_ss_hsa_mw_synergies_50sim = emba::get_synergy_scores(topo_ss_hsa_mw_50sim_file, file_type = "modelwise")
topo_ss_hsa_ew_synergies_150sim = emba::get_synergy_scores(topo_ss_hsa_ew_150sim_file)
topo_ss_hsa_mw_synergies_150sim = emba::get_synergy_scores(topo_ss_hsa_mw_150sim_file, file_type = "modelwise")

## HSA results rand
topo_rand_hsa_ew_50sim_file = paste0("results/topology-only/cascade_2.0_rand_50sim_fixpoints_hsa_ensemblewise_synergies.tab")
topo_rand_hsa_mw_50sim_file = paste0("results/topology-only/cascade_2.0_rand_50sim_fixpoints_hsa_modelwise_synergies.tab")
topo_rand_hsa_ew_150sim_file = paste0("results/topology-only/cascade_2.0_rand_150sim_fixpoints_hsa_ensemblewise_synergies.tab")
topo_rand_hsa_mw_150sim_file = paste0("results/topology-only/cascade_2.0_rand_150sim_fixpoints_hsa_modelwise_synergies.tab")

topo_rand_hsa_ew_synergies_50sim = emba::get_synergy_scores(topo_rand_hsa_ew_50sim_file)
topo_rand_hsa_mw_synergies_50sim = emba::get_synergy_scores(topo_rand_hsa_mw_50sim_file, file_type = "modelwise")
topo_rand_hsa_ew_synergies_150sim = emba::get_synergy_scores(topo_rand_hsa_ew_150sim_file)
topo_rand_hsa_mw_synergies_150sim = emba::get_synergy_scores(topo_rand_hsa_mw_150sim_file, file_type = "modelwise")

## Bliss results ss
topo_ss_bliss_ew_50sim_file = paste0("results/topology-only/cascade_2.0_ss_50sim_fixpoints_bliss_ensemblewise_synergies.tab")
topo_ss_bliss_mw_50sim_file = paste0("results/topology-only/cascade_2.0_ss_50sim_fixpoints_bliss_modelwise_synergies.tab")
topo_ss_bliss_ew_150sim_file = paste0("results/topology-only/cascade_2.0_ss_150sim_fixpoints_bliss_ensemblewise_synergies.tab")
topo_ss_bliss_mw_150sim_file = paste0("results/topology-only/cascade_2.0_ss_150sim_fixpoints_bliss_modelwise_synergies.tab")

topo_ss_bliss_ew_synergies_50sim = emba::get_synergy_scores(topo_ss_bliss_ew_50sim_file)
topo_ss_bliss_mw_synergies_50sim = emba::get_synergy_scores(topo_ss_bliss_mw_50sim_file, file_type = "modelwise")
topo_ss_bliss_ew_synergies_150sim = emba::get_synergy_scores(topo_ss_bliss_ew_150sim_file)
topo_ss_bliss_mw_synergies_150sim = emba::get_synergy_scores(topo_ss_bliss_mw_150sim_file, file_type = "modelwise")

## Bliss results rand
topo_rand_bliss_ew_50sim_file = paste0("results/topology-only/cascade_2.0_rand_50sim_fixpoints_bliss_ensemblewise_synergies.tab")
topo_rand_bliss_mw_50sim_file = paste0("results/topology-only/cascade_2.0_rand_50sim_fixpoints_bliss_modelwise_synergies.tab")
topo_rand_bliss_ew_150sim_file = paste0("results/topology-only/cascade_2.0_rand_150sim_fixpoints_bliss_ensemblewise_synergies.tab")
topo_rand_bliss_mw_150sim_file = paste0("results/topology-only/cascade_2.0_rand_150sim_fixpoints_bliss_modelwise_synergies.tab")

topo_rand_bliss_ew_synergies_50sim = emba::get_synergy_scores(topo_rand_bliss_ew_50sim_file)
topo_rand_bliss_mw_synergies_50sim = emba::get_synergy_scores(topo_rand_bliss_mw_50sim_file, file_type = "modelwise")
topo_rand_bliss_ew_synergies_150sim = emba::get_synergy_scores(topo_rand_bliss_ew_150sim_file)
topo_rand_bliss_mw_synergies_150sim = emba::get_synergy_scores(topo_rand_bliss_mw_150sim_file, file_type = "modelwise")

# calculate probability of synergy in the modelwise results
topo_ss_hsa_mw_synergies_50sim = topo_ss_hsa_mw_synergies_50sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topo_ss_hsa_mw_synergies_150sim = topo_ss_hsa_mw_synergies_150sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topo_rand_hsa_mw_synergies_50sim = topo_rand_hsa_mw_synergies_50sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topo_rand_hsa_mw_synergies_150sim = topo_rand_hsa_mw_synergies_150sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topo_ss_bliss_mw_synergies_50sim = topo_ss_bliss_mw_synergies_50sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topo_ss_bliss_mw_synergies_150sim = topo_ss_bliss_mw_synergies_150sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topo_rand_bliss_mw_synergies_50sim = topo_rand_bliss_mw_synergies_50sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topo_rand_bliss_mw_synergies_150sim = topo_rand_bliss_mw_synergies_150sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))

# Tidy the data
pred_topo_ew_hsa = bind_cols(
  topo_ss_hsa_ew_synergies_50sim %>% rename(ss_score_50sim = score),
  topo_ss_hsa_ew_synergies_150sim %>% select(score) %>% rename(ss_score_150sim = score),
  topo_rand_hsa_ew_synergies_50sim %>% select(score) %>% rename(rand_score_50sim = score),
  topo_rand_hsa_ew_synergies_150sim %>% select(score) %>% rename(rand_score_150sim = score),
  as_tibble_col(observed, column_name = "observed"))

pred_topo_mw_hsa = bind_cols(
  topo_ss_hsa_mw_synergies_50sim %>% select(perturbation, synergy_prob_ss) %>% rename(synergy_prob_ss_50sim = synergy_prob_ss),
  topo_ss_hsa_mw_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_150sim = synergy_prob_ss),
  topo_rand_hsa_mw_synergies_50sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_rand_50sim = synergy_prob_ss),
  topo_rand_hsa_mw_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_rand_150sim = synergy_prob_ss),
  as_tibble_col(observed, column_name = "observed"))

pred_topo_ew_bliss = bind_cols(
  topo_ss_bliss_ew_synergies_50sim %>% rename(ss_score_50sim = score),
  topo_ss_bliss_ew_synergies_150sim %>% select(score) %>% rename(ss_score_150sim = score),
  topo_rand_bliss_ew_synergies_50sim %>% select(score) %>% rename(rand_score_50sim = score),
  topo_rand_bliss_ew_synergies_150sim %>% select(score) %>% rename(rand_score_150sim = score),
  as_tibble_col(observed, column_name = "observed"))

pred_topo_mw_bliss = bind_cols(
  topo_ss_bliss_mw_synergies_50sim %>% select(perturbation, synergy_prob_ss) %>% rename(synergy_prob_ss_50sim = synergy_prob_ss),
  topo_ss_bliss_mw_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_150sim = synergy_prob_ss),
  topo_rand_bliss_mw_synergies_50sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_rand_50sim = synergy_prob_ss),
  topo_rand_bliss_mw_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_rand_150sim = synergy_prob_ss),
  as_tibble_col(observed, column_name = "observed"))
```

### ROC curves (HSA) {-}


```r
topo_res_ss_ew_50sim = get_roc_stats(df = pred_topo_ew_hsa, pred_col = "ss_score_50sim", label_col = "observed")
topo_res_ss_ew_150sim = get_roc_stats(df = pred_topo_ew_hsa, pred_col = "ss_score_150sim", label_col = "observed")
topo_res_rand_ew_50sim = get_roc_stats(df = pred_topo_ew_hsa, pred_col = "rand_score_50sim", label_col = "observed")
topo_res_rand_ew_150sim = get_roc_stats(df = pred_topo_ew_hsa, pred_col = "rand_score_150sim", label_col = "observed")

topo_res_ss_mw_50sim = get_roc_stats(df = pred_topo_mw_hsa, pred_col = "synergy_prob_ss_50sim", label_col = "observed", direction = ">")
topo_res_ss_mw_150sim = get_roc_stats(df = pred_topo_mw_hsa, pred_col = "synergy_prob_ss_150sim", label_col = "observed", direction = ">")
topo_res_rand_mw_50sim = get_roc_stats(df = pred_topo_mw_hsa, pred_col = "synergy_prob_rand_50sim", label_col = "observed", direction = ">")
topo_res_rand_mw_150sim = get_roc_stats(df = pred_topo_mw_hsa, pred_col = "synergy_prob_rand_150sim", label_col = "observed", direction = ">")

# Plot ROCs
plot(x = topo_res_ss_ew_50sim$roc_stats$FPR, y = topo_res_ss_ew_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (HSA)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = topo_res_ss_ew_150sim$roc_stats$FPR, y = topo_res_ss_ew_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = topo_res_rand_ew_50sim$roc_stats$FPR, y = topo_res_rand_ew_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = topo_res_rand_ew_150sim$roc_stats$FPR, y = topo_res_rand_ew_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
legend('bottomright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(topo_res_ss_ew_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(topo_res_ss_ew_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(topo_res_rand_ew_50sim$AUC, digits = 2), "Proliferative (50 sim)"),
    paste(round(topo_res_rand_ew_150sim$AUC, digits = 2), "Proliferative (150 sim)")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(x = topo_res_ss_mw_50sim$roc_stats$FPR, y = topo_res_ss_mw_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Model-wise synergies (HSA)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = topo_res_ss_mw_150sim$roc_stats$FPR, y = topo_res_ss_mw_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = topo_res_rand_mw_50sim$roc_stats$FPR, y = topo_res_rand_mw_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = topo_res_rand_mw_150sim$roc_stats$FPR, y = topo_res_rand_mw_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
legend('bottomright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(topo_res_ss_mw_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(topo_res_ss_mw_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(topo_res_rand_mw_50sim$AUC, digits = 2), "Proliferative (50 sim)"),
    paste(round(topo_res_rand_mw_150sim$AUC, digits = 2), "Proliferative (150 sim)")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
```

<img src="index_files/figure-html/ROC curves Topology Mutations (Cascade 2.0 - HSA)-1.png" width="50%" /><img src="index_files/figure-html/ROC curves Topology Mutations (Cascade 2.0 - HSA)-2.png" width="50%" />

### PR curves (HSA) {-}


```r
pr_topo_res_ss_ew_50sim = pr.curve(scores.class0 = pred_topo_ew_hsa %>% pull(ss_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topo_ew_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_topo_res_ss_ew_150sim = pr.curve(scores.class0 = pred_topo_ew_hsa %>% pull(ss_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topo_ew_hsa %>% pull(observed), curve = TRUE)
pr_topo_res_rand_ew_50sim = pr.curve(scores.class0 = pred_topo_ew_hsa %>% pull(rand_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topo_ew_hsa %>% pull(observed), curve = TRUE)
pr_topo_res_rand_ew_150sim = pr.curve(scores.class0 = pred_topo_ew_hsa %>% pull(rand_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topo_ew_hsa %>% pull(observed), curve = TRUE)

pr_topo_res_ss_mw_50sim = pr.curve(scores.class0 = pred_topo_mw_hsa %>% pull(synergy_prob_ss_50sim),
  weights.class0 = pred_topo_mw_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_topo_res_ss_mw_150sim = pr.curve(scores.class0 = pred_topo_mw_hsa %>% pull(synergy_prob_ss_150sim),
  weights.class0 = pred_topo_mw_hsa %>% pull(observed), curve = TRUE)
pr_topo_res_rand_mw_50sim = pr.curve(scores.class0 = pred_topo_mw_hsa %>% pull(synergy_prob_rand_50sim),
  weights.class0 = pred_topo_mw_hsa %>% pull(observed), curve = TRUE)
pr_topo_res_rand_mw_150sim = pr.curve(scores.class0 = pred_topo_mw_hsa %>% pull(synergy_prob_rand_150sim),
  weights.class0 = pred_topo_mw_hsa %>% pull(observed), curve = TRUE)

plot(pr_topo_res_ss_ew_50sim, main = 'PR curve, Ensemble-wise synergies (HSA)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_topo_res_ss_ew_150sim, add = TRUE, color = my_palette[2])
plot(pr_topo_res_rand_ew_50sim, add = TRUE, color = my_palette[3])
plot(pr_topo_res_rand_ew_150sim, add = TRUE, color = my_palette[4])
legend('topright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(pr_topo_res_ss_ew_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_topo_res_ss_ew_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_topo_res_rand_ew_50sim$auc.davis.goadrich, digits = 3), "Proliferative (50 sim)"),
    paste(round(pr_topo_res_rand_ew_150sim$auc.davis.goadrich, digits = 3), "Proliferative (150 sim)")))
grid(lwd = 0.5)

plot(pr_topo_res_ss_mw_50sim, main = 'PR curve, Model-wise synergies (HSA)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_topo_res_ss_mw_150sim, add = TRUE, color = my_palette[2])
plot(pr_topo_res_rand_mw_50sim, add = TRUE, color = my_palette[3])
plot(pr_topo_res_rand_mw_150sim, add = TRUE, color = my_palette[4])
legend('topright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(pr_topo_res_ss_mw_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_topo_res_ss_mw_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_topo_res_rand_mw_50sim$auc.davis.goadrich, digits = 3), "Proliferative (50 sim)"),
    paste(round(pr_topo_res_rand_mw_150sim$auc.davis.goadrich, digits = 3), "Proliferative (150 sim)")))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/PR curves Topology Mutations (Cascade 2.0 - HSA)-1.png" width="50%" /><img src="index_files/figure-html/PR curves Topology Mutations (Cascade 2.0 - HSA)-2.png" width="50%" />

### ROC-AUC sensitivity (HSA) {-}

Combine the $150$ simulation results (calibrated + proliferative)

```r
# Ensemble-wise
betas = seq(from = -20, to = 20, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_topo_ew_hsa = pred_topo_ew_hsa %>% mutate(combined_score = ss_score_150sim + beta * rand_score_150sim)
  res = get_roc_stats(df = pred_topo_ew_hsa, pred_col = "combined_score", label_col = "observed")
  auc_value = res$AUC
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2]) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/ROC-AUC sensitivity Topology Mutations (HSA - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

auc_values_mw = sapply(weights, function(w) {
  pred_topo_mw_hsa = pred_topo_mw_hsa %>% 
    mutate(weighted_prob = (1 - w) * pred_topo_mw_hsa$synergy_prob_ss_150sim + w * pred_topo_mw_hsa$synergy_prob_rand_150sim)
  res = get_roc_stats(df = pred_topo_mw_hsa, pred_col = "weighted_prob", label_col = "observed", direction = ">")
  auc_value = res$AUC
})

df_mw = as_tibble(cbind(weights, auc_values_mw))

ggline(data = df_mw, x = "weights", y = "auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/ROC-AUC sensitivity Topology Mutations (HSA - Cascade 2.0)-2.png" width="80%" style="display: block; margin: auto;" />

### PR-AUC sensitivity (HSA) {-}

Combine the $150$ simulation results (calibrated + proliferative)

```r
# Ensemble-wise
betas = seq(from = -10, to = 10, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_topo_ew_hsa = pred_topo_ew_hsa %>% mutate(combined_score = ss_score_150sim + beta * rand_score_150sim)
  res = pr.curve(scores.class0 = pred_topo_ew_hsa %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_topo_ew_hsa %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under PR Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2], label = "betas", repel = TRUE,
  label.select = list(criteria = "`y` == max(`y`)"),
  font.label = list(size = 14)) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/PR-AUC sensitivity Topology Mutations (HSA - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

### ROC curves (Bliss) {-}


```r
topo_res_ss_ew_50sim = get_roc_stats(df = pred_topo_ew_bliss, pred_col = "ss_score_50sim", label_col = "observed")
topo_res_ss_ew_150sim = get_roc_stats(df = pred_topo_ew_bliss, pred_col = "ss_score_150sim", label_col = "observed")
topo_res_rand_ew_50sim = get_roc_stats(df = pred_topo_ew_bliss, pred_col = "rand_score_50sim", label_col = "observed")
topo_res_rand_ew_150sim = get_roc_stats(df = pred_topo_ew_bliss, pred_col = "rand_score_150sim", label_col = "observed")

topo_res_ss_mw_50sim = get_roc_stats(df = pred_topo_mw_bliss, pred_col = "synergy_prob_ss_50sim", label_col = "observed", direction = ">")
topo_res_ss_mw_150sim = get_roc_stats(df = pred_topo_mw_bliss, pred_col = "synergy_prob_ss_150sim", label_col = "observed", direction = ">")
topo_res_rand_mw_50sim = get_roc_stats(df = pred_topo_mw_bliss, pred_col = "synergy_prob_rand_50sim", label_col = "observed", direction = ">")
topo_res_rand_mw_150sim = get_roc_stats(df = pred_topo_mw_bliss, pred_col = "synergy_prob_rand_150sim", label_col = "observed", direction = ">")

# Plot ROCs
plot(x = topo_res_ss_ew_50sim$roc_stats$FPR, y = topo_res_ss_ew_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = topo_res_ss_ew_150sim$roc_stats$FPR, y = topo_res_ss_ew_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = topo_res_rand_ew_50sim$roc_stats$FPR, y = topo_res_rand_ew_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = topo_res_rand_ew_150sim$roc_stats$FPR, y = topo_res_rand_ew_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
legend('bottomright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(topo_res_ss_ew_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(topo_res_ss_ew_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(topo_res_rand_ew_50sim$AUC, digits = 2), "Proliferative (50 sim)"),
    paste(round(topo_res_rand_ew_150sim$AUC, digits = 2), "Proliferative (150 sim)")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(x = topo_res_ss_mw_50sim$roc_stats$FPR, y = topo_res_ss_mw_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Model-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = topo_res_ss_mw_150sim$roc_stats$FPR, y = topo_res_ss_mw_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = topo_res_rand_mw_50sim$roc_stats$FPR, y = topo_res_rand_mw_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = topo_res_rand_mw_150sim$roc_stats$FPR, y = topo_res_rand_mw_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
legend('bottomright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(topo_res_ss_mw_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(topo_res_ss_mw_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(topo_res_rand_mw_50sim$AUC, digits = 2), "Proliferative (50 sim)"),
    paste(round(topo_res_rand_mw_150sim$AUC, digits = 2), "Proliferative (150 sim)")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
```

<img src="index_files/figure-html/ROC curves Topology Mutations (Cascade 2.0 - Bliss)-1.png" width="50%" /><img src="index_files/figure-html/ROC curves Topology Mutations (Cascade 2.0 - Bliss)-2.png" width="50%" />

### PR curves (Bliss) {-}


```r
pr_topo_res_ss_ew_50sim = pr.curve(scores.class0 = pred_topo_ew_bliss %>% pull(ss_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topo_ew_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_topo_res_ss_ew_150sim = pr.curve(scores.class0 = pred_topo_ew_bliss %>% pull(ss_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topo_ew_bliss %>% pull(observed), curve = TRUE)
pr_topo_res_rand_ew_50sim = pr.curve(scores.class0 = pred_topo_ew_bliss %>% pull(rand_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topo_ew_bliss %>% pull(observed), curve = TRUE)
pr_topo_res_rand_ew_150sim = pr.curve(scores.class0 = pred_topo_ew_bliss %>% pull(rand_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topo_ew_bliss %>% pull(observed), curve = TRUE)

pr_topo_res_ss_mw_50sim = pr.curve(scores.class0 = pred_topo_mw_bliss %>% pull(synergy_prob_ss_50sim),
  weights.class0 = pred_topo_mw_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_topo_res_ss_mw_150sim = pr.curve(scores.class0 = pred_topo_mw_bliss %>% pull(synergy_prob_ss_150sim),
  weights.class0 = pred_topo_mw_bliss %>% pull(observed), curve = TRUE)
pr_topo_res_rand_mw_50sim = pr.curve(scores.class0 = pred_topo_mw_bliss %>% pull(synergy_prob_rand_50sim),
  weights.class0 = pred_topo_mw_bliss %>% pull(observed), curve = TRUE)
pr_topo_res_rand_mw_150sim = pr.curve(scores.class0 = pred_topo_mw_bliss %>% pull(synergy_prob_rand_150sim),
  weights.class0 = pred_topo_mw_bliss %>% pull(observed), curve = TRUE)

plot(pr_topo_res_ss_ew_50sim, main = 'PR curve, Ensemble-wise synergies (Bliss)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_topo_res_ss_ew_150sim, add = TRUE, color = my_palette[2])
plot(pr_topo_res_rand_ew_50sim, add = TRUE, color = my_palette[3])
plot(pr_topo_res_rand_ew_150sim, add = TRUE, color = my_palette[4])
legend('topright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(pr_topo_res_ss_ew_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_topo_res_ss_ew_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_topo_res_rand_ew_50sim$auc.davis.goadrich, digits = 3), "Proliferative (50 sim)"),
    paste(round(pr_topo_res_rand_ew_150sim$auc.davis.goadrich, digits = 3), "Proliferative (150 sim)")))
grid(lwd = 0.5)

plot(pr_topo_res_ss_mw_50sim, main = 'PR curve, Model-wise synergies (Bliss)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_topo_res_ss_mw_150sim, add = TRUE, color = my_palette[2])
plot(pr_topo_res_rand_mw_50sim, add = TRUE, color = my_palette[3])
plot(pr_topo_res_rand_mw_150sim, add = TRUE, color = my_palette[4])
legend('topright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(pr_topo_res_ss_mw_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_topo_res_ss_mw_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_topo_res_rand_mw_50sim$auc.davis.goadrich, digits = 3), "Proliferative (50 sim)"),
    paste(round(pr_topo_res_rand_mw_150sim$auc.davis.goadrich, digits = 3), "Proliferative (150 sim)")))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/PR curves Topology Mutations (Cascade 2.0 - Bliss)-1.png" width="50%" /><img src="index_files/figure-html/PR curves Topology Mutations (Cascade 2.0 - Bliss)-2.png" width="50%" />

### ROC-AUC sensitivity (Bliss) {-}


```r
betas = seq(from = -10, to = 10, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_topo_ew_bliss = pred_topo_ew_bliss %>% mutate(combined_score = ss_score_150sim + beta * rand_score_150sim)
  res = get_roc_stats(df = pred_topo_ew_bliss, pred_col = "combined_score", label_col = "observed")
  auc_value = res$AUC
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2], label = "betas", repel = TRUE,
  label.select = list(criteria = "`y` == max(`y`)"),
  font.label = list(size = 14)) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/ROC-AUC sensitivity Topology Mutations (Bliss - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

### PR-AUC sensitivity (Bliss) {-}


```r
# Ensemble-wise
betas = seq(from = -10, to = 10, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_topo_ew_bliss = pred_topo_ew_bliss %>% mutate(combined_score = ss_score_150sim + beta * rand_score_150sim)
  res = pr.curve(scores.class0 = pred_topo_ew_bliss %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_topo_ew_bliss %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under PR Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2], label = "betas", repel = TRUE,
  label.select = list(criteria = "`y` == max(`y`)"),
  font.label = list(size = 14)) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/PR-AUC sensitivity Topology Mutations (Bliss - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

### ROC and PRC for best beta (Bliss) {-}


```r
best_beta = -1
pred_topo_ew_bliss = pred_topo_ew_bliss %>% mutate(best_score = ss_score_150sim + best_beta * rand_score_150sim)

roc_best_res = get_roc_stats(df = pred_topo_ew_bliss, pred_col = "best_score", label_col = "observed")
pr_best_res = pr.curve(scores.class0 = pred_topo_ew_bliss %>% pull(best_score) %>% (function(x) {-x}), 
    weights.class0 = pred_topo_ew_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)

# Plot best ROC
plot(x = roc_best_res$roc_stats$FPR, y = roc_best_res$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = TeX('ROC curve (Ensemble-wise), $calibrated + \\beta \\times proliferative$'),
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
legend('bottomright', title = TeX('AUC ($\\beta$ = -1)'), col = my_palette[1], pch = 19,
  legend = paste(round(roc_best_res$AUC, digits = 2), 'Bliss (150 sim)'))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

# Plot best PRC
plot(pr_best_res, main = TeX('PR curve (Ensemble-wise), $calibrated + \\beta \\times proliferative$'),
  auc.main = FALSE, color = my_palette[2], rand.plot = TRUE)
legend('topright', title = TeX('AUC ($\\beta$ = -1)'), col = my_palette[2], pch = 19,
  legend = paste(round(pr_best_res$auc.davis.goadrich, digits = 3), 'Bliss (150 sim)'))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/ROC and PRC for best beta - Topology Mutations (Bliss - Cascade 2.0)-1.png" width="50%" /><img src="index_files/figure-html/ROC and PRC for best beta - Topology Mutations (Bliss - Cascade 2.0)-2.png" width="50%" />

## Topology and Link Operator Mutations {-}

:::{.note}
We run `Gitsbe` simulations with $50$ topology mutations and $3000$ link operator mutations (bootstrap values, reduced to $10$ and $3$ respectively after models with stabla states have been found), both for $50$ and $150$ simulations and both **fitting to steady state** (calibrated models) and to a **proliferative phenotype** (so not random models but as close as it can get to that since we are discussing topology mutations).

Results were assessed via `Drabme`'s ensemble and modelwise results, using both the HSA and Bliss method.
:::


```r
# 'ss' => calibrated models, 'rand' => proliferative models (so not random but kind of!)
# 'ew' => ensemble-wise, 'mw' => modelwise

## HSA results ss
topolink_ss_hsa_ew_50sim_file = paste0("results/topo-and-link/cascade_2.0_ss_50sim_fixpoints_hsa_ensemblewise_synergies.tab")
topolink_ss_hsa_mw_50sim_file = paste0("results/topo-and-link/cascade_2.0_ss_50sim_fixpoints_hsa_modelwise_synergies.tab")
topolink_ss_hsa_ew_150sim_file = paste0("results/topo-and-link/cascade_2.0_ss_150sim_fixpoints_hsa_ensemblewise_synergies.tab")
topolink_ss_hsa_mw_150sim_file = paste0("results/topo-and-link/cascade_2.0_ss_150sim_fixpoints_hsa_modelwise_synergies.tab")

topolink_ss_hsa_ew_synergies_50sim = emba::get_synergy_scores(topolink_ss_hsa_ew_50sim_file)
topolink_ss_hsa_mw_synergies_50sim = emba::get_synergy_scores(topolink_ss_hsa_mw_50sim_file, file_type = "modelwise")
topolink_ss_hsa_ew_synergies_150sim = emba::get_synergy_scores(topolink_ss_hsa_ew_150sim_file)
topolink_ss_hsa_mw_synergies_150sim = emba::get_synergy_scores(topolink_ss_hsa_mw_150sim_file, file_type = "modelwise")

## HSA results rand
topolink_rand_hsa_ew_50sim_file = paste0("results/topo-and-link/cascade_2.0_rand_50sim_fixpoints_hsa_ensemblewise_synergies.tab")
topolink_rand_hsa_mw_50sim_file = paste0("results/topo-and-link/cascade_2.0_rand_50sim_fixpoints_hsa_modelwise_synergies.tab")
topolink_rand_hsa_ew_150sim_file = paste0("results/topo-and-link/cascade_2.0_rand_150sim_fixpoints_hsa_ensemblewise_synergies.tab")
topolink_rand_hsa_mw_150sim_file = paste0("results/topo-and-link/cascade_2.0_rand_150sim_fixpoints_hsa_modelwise_synergies.tab")

topolink_rand_hsa_ew_synergies_50sim = emba::get_synergy_scores(topolink_rand_hsa_ew_50sim_file)
topolink_rand_hsa_mw_synergies_50sim = emba::get_synergy_scores(topolink_rand_hsa_mw_50sim_file, file_type = "modelwise")
topolink_rand_hsa_ew_synergies_150sim = emba::get_synergy_scores(topolink_rand_hsa_ew_150sim_file)
topolink_rand_hsa_mw_synergies_150sim = emba::get_synergy_scores(topolink_rand_hsa_mw_150sim_file, file_type = "modelwise")

## Bliss results ss
topolink_ss_bliss_ew_50sim_file = paste0("results/topo-and-link/cascade_2.0_ss_50sim_fixpoints_bliss_ensemblewise_synergies.tab")
topolink_ss_bliss_mw_50sim_file = paste0("results/topo-and-link/cascade_2.0_ss_50sim_fixpoints_bliss_modelwise_synergies.tab")
topolink_ss_bliss_ew_150sim_file = paste0("results/topo-and-link/cascade_2.0_ss_150sim_fixpoints_bliss_ensemblewise_synergies.tab")
topolink_ss_bliss_mw_150sim_file = paste0("results/topo-and-link/cascade_2.0_ss_150sim_fixpoints_bliss_modelwise_synergies.tab")

topolink_ss_bliss_ew_synergies_50sim = emba::get_synergy_scores(topolink_ss_bliss_ew_50sim_file)
topolink_ss_bliss_mw_synergies_50sim = emba::get_synergy_scores(topolink_ss_bliss_mw_50sim_file, file_type = "modelwise")
topolink_ss_bliss_ew_synergies_150sim = emba::get_synergy_scores(topolink_ss_bliss_ew_150sim_file)
topolink_ss_bliss_mw_synergies_150sim = emba::get_synergy_scores(topolink_ss_bliss_mw_150sim_file, file_type = "modelwise")

## Bliss results rand
topolink_rand_bliss_ew_50sim_file = paste0("results/topo-and-link/cascade_2.0_rand_50sim_fixpoints_bliss_ensemblewise_synergies.tab")
topolink_rand_bliss_mw_50sim_file = paste0("results/topo-and-link/cascade_2.0_rand_50sim_fixpoints_bliss_modelwise_synergies.tab")
topolink_rand_bliss_ew_150sim_file = paste0("results/topo-and-link/cascade_2.0_rand_150sim_fixpoints_bliss_ensemblewise_synergies.tab")
topolink_rand_bliss_mw_150sim_file = paste0("results/topo-and-link/cascade_2.0_rand_150sim_fixpoints_bliss_modelwise_synergies.tab")

topolink_rand_bliss_ew_synergies_50sim = emba::get_synergy_scores(topolink_rand_bliss_ew_50sim_file)
topolink_rand_bliss_mw_synergies_50sim = emba::get_synergy_scores(topolink_rand_bliss_mw_50sim_file, file_type = "modelwise")
topolink_rand_bliss_ew_synergies_150sim = emba::get_synergy_scores(topolink_rand_bliss_ew_150sim_file)
topolink_rand_bliss_mw_synergies_150sim = emba::get_synergy_scores(topolink_rand_bliss_mw_150sim_file, file_type = "modelwise")

# calculate probability of synergy in the modelwise results
topolink_ss_hsa_mw_synergies_50sim = topolink_ss_hsa_mw_synergies_50sim %>% 
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topolink_ss_hsa_mw_synergies_150sim = topolink_ss_hsa_mw_synergies_150sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topolink_rand_hsa_mw_synergies_50sim = topolink_rand_hsa_mw_synergies_50sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topolink_rand_hsa_mw_synergies_150sim = topolink_rand_hsa_mw_synergies_150sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topolink_ss_bliss_mw_synergies_50sim = topolink_ss_bliss_mw_synergies_50sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topolink_ss_bliss_mw_synergies_150sim = topolink_ss_bliss_mw_synergies_150sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topolink_rand_bliss_mw_synergies_50sim = topolink_rand_bliss_mw_synergies_50sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))
topolink_rand_bliss_mw_synergies_150sim = topolink_rand_bliss_mw_synergies_150sim %>%
  mutate(synergy_prob_ss = synergies/(synergies + `non-synergies`))

# Tidy the data
pred_topolink_ew_hsa = bind_cols(
  topolink_ss_hsa_ew_synergies_50sim %>% rename(ss_score_50sim = score),
  topolink_ss_hsa_ew_synergies_150sim %>% select(score) %>% rename(ss_score_150sim = score),
  topolink_rand_hsa_ew_synergies_50sim %>% select(score) %>% rename(rand_score_50sim = score),
  topolink_rand_hsa_ew_synergies_150sim %>% select(score) %>% rename(rand_score_150sim = score),
  as_tibble_col(observed, column_name = "observed"))

pred_topolink_mw_hsa = bind_cols(
  topolink_ss_hsa_mw_synergies_50sim %>% select(perturbation, synergy_prob_ss) %>% rename(synergy_prob_ss_50sim = synergy_prob_ss),
  topolink_ss_hsa_mw_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_150sim = synergy_prob_ss),
  topolink_rand_hsa_mw_synergies_50sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_rand_50sim = synergy_prob_ss),
  topolink_rand_hsa_mw_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_rand_150sim = synergy_prob_ss),
  as_tibble_col(observed, column_name = "observed"))

pred_topolink_ew_bliss = bind_cols(
  topolink_ss_bliss_ew_synergies_50sim %>% rename(ss_score_50sim = score),
  topolink_ss_bliss_ew_synergies_150sim %>% select(score) %>% rename(ss_score_150sim = score),
  topolink_rand_bliss_ew_synergies_50sim %>% select(score) %>% rename(rand_score_50sim = score),
  topolink_rand_bliss_ew_synergies_150sim %>% select(score) %>% rename(rand_score_150sim = score),
  as_tibble_col(observed, column_name = "observed"))

pred_topolink_mw_bliss = bind_cols(
  topolink_ss_bliss_mw_synergies_50sim %>% select(perturbation, synergy_prob_ss) %>% rename(synergy_prob_ss_50sim = synergy_prob_ss),
  topolink_ss_bliss_mw_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_ss_150sim = synergy_prob_ss),
  topolink_rand_bliss_mw_synergies_50sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_rand_50sim = synergy_prob_ss),
  topolink_rand_bliss_mw_synergies_150sim %>% select(synergy_prob_ss) %>% rename(synergy_prob_rand_150sim = synergy_prob_ss),
  as_tibble_col(observed, column_name = "observed"))
```

### ROC curves (HSA) {-}


```r
topolink_res_ss_ew_50sim = get_roc_stats(df = pred_topolink_ew_hsa, pred_col = "ss_score_50sim", label_col = "observed")
topolink_res_ss_ew_150sim = get_roc_stats(df = pred_topolink_ew_hsa, pred_col = "ss_score_150sim", label_col = "observed")
topolink_res_rand_ew_50sim = get_roc_stats(df = pred_topolink_ew_hsa, pred_col = "rand_score_50sim", label_col = "observed")
topolink_res_rand_ew_150sim = get_roc_stats(df = pred_topolink_ew_hsa, pred_col = "rand_score_150sim", label_col = "observed")

topolink_res_ss_mw_50sim = get_roc_stats(df = pred_topolink_mw_hsa, pred_col = "synergy_prob_ss_50sim", label_col = "observed", direction = ">")
topolink_res_ss_mw_150sim = get_roc_stats(df = pred_topolink_mw_hsa, pred_col = "synergy_prob_ss_150sim", label_col = "observed", direction = ">")
topolink_res_rand_mw_50sim = get_roc_stats(df = pred_topolink_mw_hsa, pred_col = "synergy_prob_rand_50sim", label_col = "observed", direction = ">")
topolink_res_rand_mw_150sim = get_roc_stats(df = pred_topolink_mw_hsa, pred_col = "synergy_prob_rand_150sim", label_col = "observed", direction = ">")

# Plot ROCs
plot(x = topolink_res_ss_ew_50sim$roc_stats$FPR, y = topolink_res_ss_ew_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (HSA)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = topolink_res_ss_ew_150sim$roc_stats$FPR, y = topolink_res_ss_ew_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = topolink_res_rand_ew_50sim$roc_stats$FPR, y = topolink_res_rand_ew_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = topolink_res_rand_ew_150sim$roc_stats$FPR, y = topolink_res_rand_ew_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
legend('bottomright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(topolink_res_ss_ew_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(topolink_res_ss_ew_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(topolink_res_rand_ew_50sim$AUC, digits = 2), "Proliferative (50 sim)"),
    paste(round(topolink_res_rand_ew_150sim$AUC, digits = 2), "Proliferative (150 sim)")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(x = topolink_res_ss_mw_50sim$roc_stats$FPR, y = topolink_res_ss_mw_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Model-wise synergies (HSA)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = topolink_res_ss_mw_150sim$roc_stats$FPR, y = topolink_res_ss_mw_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = topolink_res_rand_mw_50sim$roc_stats$FPR, y = topolink_res_rand_mw_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = topolink_res_rand_mw_150sim$roc_stats$FPR, y = topolink_res_rand_mw_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
legend('bottomright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(topolink_res_ss_mw_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(topolink_res_ss_mw_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(topolink_res_rand_mw_50sim$AUC, digits = 2), "Proliferative (50 sim)"),
    paste(round(topolink_res_rand_mw_150sim$AUC, digits = 2), "Proliferative (150 sim)")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
```

<img src="index_files/figure-html/ROC curves Topology and Link Operator Mutations (Cascade 2.0 - HSA)-1.png" width="50%" /><img src="index_files/figure-html/ROC curves Topology and Link Operator Mutations (Cascade 2.0 - HSA)-2.png" width="50%" />

### PR curves (HSA) {-}


```r
pr_topolink_res_ss_ew_50sim = pr.curve(scores.class0 = pred_topolink_ew_hsa %>% pull(ss_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topolink_ew_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_topolink_res_ss_ew_150sim = pr.curve(scores.class0 = pred_topolink_ew_hsa %>% pull(ss_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topolink_ew_hsa %>% pull(observed), curve = TRUE)
pr_topolink_res_rand_ew_50sim = pr.curve(scores.class0 = pred_topolink_ew_hsa %>% pull(rand_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topolink_ew_hsa %>% pull(observed), curve = TRUE)
pr_topolink_res_rand_ew_150sim = pr.curve(scores.class0 = pred_topolink_ew_hsa %>% pull(rand_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topolink_ew_hsa %>% pull(observed), curve = TRUE)

pr_topolink_res_ss_mw_50sim = pr.curve(scores.class0 = pred_topolink_mw_hsa %>% pull(synergy_prob_ss_50sim),
  weights.class0 = pred_topolink_mw_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_topolink_res_ss_mw_150sim = pr.curve(scores.class0 = pred_topolink_mw_hsa %>% pull(synergy_prob_ss_150sim),
  weights.class0 = pred_topolink_mw_hsa %>% pull(observed), curve = TRUE)
pr_topolink_res_rand_mw_50sim = pr.curve(scores.class0 = pred_topolink_mw_hsa %>% pull(synergy_prob_rand_50sim),
  weights.class0 = pred_topolink_mw_hsa %>% pull(observed), curve = TRUE)
pr_topolink_res_rand_mw_150sim = pr.curve(scores.class0 = pred_topolink_mw_hsa %>% pull(synergy_prob_rand_150sim),
  weights.class0 = pred_topolink_mw_hsa %>% pull(observed), curve = TRUE)

plot(pr_topolink_res_ss_ew_50sim, main = 'PR curve, Ensemble-wise synergies (HSA)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_topolink_res_ss_ew_150sim, add = TRUE, color = my_palette[2])
plot(pr_topolink_res_rand_ew_50sim, add = TRUE, color = my_palette[3])
plot(pr_topolink_res_rand_ew_150sim, add = TRUE, color = my_palette[4])
legend('topright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(pr_topolink_res_ss_ew_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_topolink_res_ss_ew_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_topolink_res_rand_ew_50sim$auc.davis.goadrich, digits = 3), "Proliferative (50 sim)"),
    paste(round(pr_topolink_res_rand_ew_150sim$auc.davis.goadrich, digits = 3), "Proliferative (150 sim)")))
grid(lwd = 0.5)

plot(pr_topolink_res_ss_mw_50sim, main = 'PR curve, Model-wise synergies (HSA)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_topolink_res_ss_mw_150sim, add = TRUE, color = my_palette[2])
plot(pr_topolink_res_rand_mw_50sim, add = TRUE, color = my_palette[3])
plot(pr_topolink_res_rand_mw_150sim, add = TRUE, color = my_palette[4])
legend('topright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(pr_topolink_res_ss_mw_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_topolink_res_ss_mw_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_topolink_res_rand_mw_50sim$auc.davis.goadrich, digits = 3), "Proliferative (50 sim)"),
    paste(round(pr_topolink_res_rand_mw_150sim$auc.davis.goadrich, digits = 3), "Proliferative (150 sim)")))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/PR curves Topology and Link Mutations (Cascade 2.0 - HSA)-1.png" width="50%" /><img src="index_files/figure-html/PR curves Topology and Link Mutations (Cascade 2.0 - HSA)-2.png" width="50%" />

### ROC-AUC sensitivity (HSA) {-}

Combine the $150$ simulation results (calibrated + proliferative)

```r
# Ensemble-wise
betas = seq(from = -10, to = 10, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_topolink_ew_hsa = pred_topolink_ew_hsa %>% mutate(combined_score = ss_score_150sim + beta * rand_score_150sim)
  res = get_roc_stats(df = pred_topolink_ew_hsa, pred_col = "combined_score", label_col = "observed")
  auc_value = res$AUC
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2], label = "betas", label.rectangle = TRUE, repel = TRUE,
  label.select = list(criteria = "`y` == max(`y`)"),
  font.label = list(size = 14)) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/ROC-AUC sensitivity Topology and Link Mutations (HSA - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

auc_values_mw = sapply(weights, function(w) {
  pred_topolink_mw_hsa = pred_topolink_mw_hsa %>% 
    mutate(weighted_prob = (1 - w) * pred_topolink_mw_hsa$synergy_prob_ss_150sim + w * pred_topolink_mw_hsa$synergy_prob_rand_150sim)
  res = get_roc_stats(df = pred_topolink_mw_hsa, pred_col = "weighted_prob", label_col = "observed", direction = ">")
  auc_value = res$AUC
})

df_mw = as_tibble(cbind(weights, auc_values_mw))

ggline(data = df_mw, x = "weights", y = "auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/ROC-AUC sensitivity Topology and Link Mutations (HSA - Cascade 2.0)-2.png" width="80%" style="display: block; margin: auto;" />

### PR-AUC sensitivity (HSA) {-}

Combine the $150$ simulation results (calibrated + proliferative)

```r
# Ensemble-wise
betas = seq(from = -10, to = 10, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_topolink_ew_hsa = pred_topolink_ew_hsa %>% mutate(combined_score = ss_score_150sim + beta * rand_score_150sim)
  res = pr.curve(scores.class0 = pred_topolink_ew_hsa %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_topolink_ew_hsa %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under PR Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2], label = "betas", label.rectangle = TRUE, repel = TRUE,
  label.select = list(criteria = "`y` == max(`y`)"),
  font.label = list(size = 14)) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/PR-AUC sensitivity Topology and Link Mutations (HSA - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

auc_values_mw = sapply(weights, function(w) {
  pred_topolink_mw_hsa = pred_topolink_mw_hsa %>% 
    mutate(weighted_prob = (1 - w) * pred_topolink_mw_hsa$synergy_prob_ss_150sim + w * pred_topolink_mw_hsa$synergy_prob_rand_150sim)
  res = pr.curve(scores.class0 = pred_topolink_mw_hsa %>% pull(weighted_prob), 
    weights.class0 = pred_topolink_mw_hsa %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_mw = as_tibble(cbind(weights, auc_values_mw))

ggline(data = df_mw, x = "weights", y = "auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "AUC (Area Under PR Curve)",
  title = TeX("AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/PR-AUC sensitivity Topology and Link Mutations (HSA - Cascade 2.0)-2.png" width="80%" style="display: block; margin: auto;" />

### ROC curves (Bliss) {-}


```r
topolink_res_ss_ew_50sim = get_roc_stats(df = pred_topolink_ew_bliss, pred_col = "ss_score_50sim", label_col = "observed")
topolink_res_ss_ew_150sim = get_roc_stats(df = pred_topolink_ew_bliss, pred_col = "ss_score_150sim", label_col = "observed")
topolink_res_rand_ew_50sim = get_roc_stats(df = pred_topolink_ew_bliss, pred_col = "rand_score_50sim", label_col = "observed")
topolink_res_rand_ew_150sim = get_roc_stats(df = pred_topolink_ew_bliss, pred_col = "rand_score_150sim", label_col = "observed")

topolink_res_ss_mw_50sim = get_roc_stats(df = pred_topolink_mw_bliss, pred_col = "synergy_prob_ss_50sim", label_col = "observed", direction = ">")
topolink_res_ss_mw_150sim = get_roc_stats(df = pred_topolink_mw_bliss, pred_col = "synergy_prob_ss_150sim", label_col = "observed", direction = ">")
topolink_res_rand_mw_50sim = get_roc_stats(df = pred_topolink_mw_bliss, pred_col = "synergy_prob_rand_50sim", label_col = "observed", direction = ">")
topolink_res_rand_mw_150sim = get_roc_stats(df = pred_topolink_mw_bliss, pred_col = "synergy_prob_rand_150sim", label_col = "observed", direction = ">")

# Plot ROCs
plot(x = topolink_res_ss_ew_50sim$roc_stats$FPR, y = topolink_res_ss_ew_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Ensemble-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = topolink_res_ss_ew_150sim$roc_stats$FPR, y = topolink_res_ss_ew_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = topolink_res_rand_ew_50sim$roc_stats$FPR, y = topolink_res_rand_ew_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = topolink_res_rand_ew_150sim$roc_stats$FPR, y = topolink_res_rand_ew_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
legend('bottomright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(topolink_res_ss_ew_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(topolink_res_ss_ew_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(topolink_res_rand_ew_50sim$AUC, digits = 2), "Proliferative (50 sim)"),
    paste(round(topolink_res_rand_ew_150sim$AUC, digits = 2), "Proliferative (150 sim)")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

plot(x = topolink_res_ss_mw_50sim$roc_stats$FPR, y = topolink_res_ss_mw_50sim$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = 'ROC curve, Model-wise synergies (Bliss)',
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = topolink_res_ss_mw_150sim$roc_stats$FPR, y = topolink_res_ss_mw_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[2])
lines(x = topolink_res_rand_mw_50sim$roc_stats$FPR, y = topolink_res_rand_mw_50sim$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
lines(x = topolink_res_rand_mw_150sim$roc_stats$FPR, y = topolink_res_rand_mw_150sim$roc_stats$TPR,
  lwd = 3, col = my_palette[4])
legend('bottomright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(topolink_res_ss_mw_50sim$AUC, digits = 2), "Calibrated (50 sim)"),
    paste(round(topolink_res_ss_mw_150sim$AUC, digits = 2), "Calibrated (150 sim)"),
    paste(round(topolink_res_rand_mw_50sim$AUC, digits = 2), "Proliferative (50 sim)"),
    paste(round(topolink_res_rand_mw_150sim$AUC, digits = 2), "Proliferative (150 sim)")))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)
```

<img src="index_files/figure-html/ROC curves Topology and Link Mutations (Cascade 2.0 - Bliss)-1.png" width="50%" /><img src="index_files/figure-html/ROC curves Topology and Link Mutations (Cascade 2.0 - Bliss)-2.png" width="50%" />

### PR curves (Bliss) {-}


```r
pr_topolink_res_ss_ew_50sim = pr.curve(scores.class0 = pred_topolink_ew_bliss %>% pull(ss_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topolink_ew_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_topolink_res_ss_ew_150sim = pr.curve(scores.class0 = pred_topolink_ew_bliss %>% pull(ss_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topolink_ew_bliss %>% pull(observed), curve = TRUE)
pr_topolink_res_rand_ew_50sim = pr.curve(scores.class0 = pred_topolink_ew_bliss %>% pull(rand_score_50sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topolink_ew_bliss %>% pull(observed), curve = TRUE)
pr_topolink_res_rand_ew_150sim = pr.curve(scores.class0 = pred_topolink_ew_bliss %>% pull(rand_score_150sim) %>% (function(x) {-x}), 
  weights.class0 = pred_topolink_ew_bliss %>% pull(observed), curve = TRUE)

pr_topolink_res_ss_mw_50sim = pr.curve(scores.class0 = pred_topolink_mw_bliss %>% pull(synergy_prob_ss_50sim),
  weights.class0 = pred_topolink_mw_bliss %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_topolink_res_ss_mw_150sim = pr.curve(scores.class0 = pred_topolink_mw_bliss %>% pull(synergy_prob_ss_150sim),
  weights.class0 = pred_topolink_mw_bliss %>% pull(observed), curve = TRUE)
pr_topolink_res_rand_mw_50sim = pr.curve(scores.class0 = pred_topolink_mw_bliss %>% pull(synergy_prob_rand_50sim),
  weights.class0 = pred_topolink_mw_bliss %>% pull(observed), curve = TRUE)
pr_topolink_res_rand_mw_150sim = pr.curve(scores.class0 = pred_topolink_mw_bliss %>% pull(synergy_prob_rand_150sim),
  weights.class0 = pred_topolink_mw_bliss %>% pull(observed), curve = TRUE)

plot(pr_topolink_res_ss_ew_50sim, main = 'PR curve, Ensemble-wise synergies (Bliss)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_topolink_res_ss_ew_150sim, add = TRUE, color = my_palette[2])
plot(pr_topolink_res_rand_ew_50sim, add = TRUE, color = my_palette[3])
plot(pr_topolink_res_rand_ew_150sim, add = TRUE, color = my_palette[4])
legend('topright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(pr_topolink_res_ss_ew_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_topolink_res_ss_ew_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_topolink_res_rand_ew_50sim$auc.davis.goadrich, digits = 3), "Proliferative (50 sim)"),
    paste(round(pr_topolink_res_rand_ew_150sim$auc.davis.goadrich, digits = 3), "Proliferative (150 sim)")))
grid(lwd = 0.5)

plot(pr_topolink_res_ss_mw_50sim, main = 'PR curve, Model-wise synergies (Bliss)',
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_topolink_res_ss_mw_150sim, add = TRUE, color = my_palette[2])
plot(pr_topolink_res_rand_mw_50sim, add = TRUE, color = my_palette[3])
plot(pr_topolink_res_rand_mw_150sim, add = TRUE, color = my_palette[4])
legend('topright', title = 'AUC', col = my_palette[1:4], pch = 19,
  legend = c(paste(round(pr_topolink_res_ss_mw_50sim$auc.davis.goadrich, digits = 3), "Calibrated (50 sim)"),
    paste(round(pr_topolink_res_ss_mw_150sim$auc.davis.goadrich, digits = 3), "Calibrated (150 sim)"),
    paste(round(pr_topolink_res_rand_mw_50sim$auc.davis.goadrich, digits = 3), "Proliferative (50 sim)"),
    paste(round(pr_topolink_res_rand_mw_150sim$auc.davis.goadrich, digits = 3), "Proliferative (150 sim)")))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/PR curves Topology and Link Mutations (Cascade 2.0 - Bliss)-1.png" width="50%" /><img src="index_files/figure-html/PR curves Topology and Link Mutations (Cascade 2.0 - Bliss)-2.png" width="50%" />

### ROC-AUC sensitivity (Bliss) {-}

Combine the $150$ simulation results (calibrated + proliferative)

```r
# Ensemble-wise
betas = seq(from = -10, to = 10, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_topolink_ew_bliss = pred_topolink_ew_bliss %>% mutate(combined_score = ss_score_150sim + beta * rand_score_150sim)
  res = get_roc_stats(df = pred_topolink_ew_bliss, pred_col = "combined_score", label_col = "observed")
  auc_value = res$AUC
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2], label = "betas", repel = TRUE,
  label.select = list(criteria = "`y` == max(`y`)"),
  font.label = list(size = 14)) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/ROC-AUC sensitivity Topology and Link Mutations (Bliss - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

auc_values_mw = sapply(weights, function(w) {
  pred_topolink_mw_bliss = pred_topolink_mw_bliss %>% 
    mutate(weighted_prob = (1 - w) * pred_topolink_mw_bliss$synergy_prob_ss_150sim + w * pred_topolink_mw_bliss$synergy_prob_rand_150sim)
  res = get_roc_stats(df = pred_topolink_mw_bliss, pred_col = "weighted_prob", label_col = "observed", direction = ">")
  auc_value = res$AUC
})

df_mw = as_tibble(cbind(weights, auc_values_mw))

ggline(data = df_mw, x = "weights", y = "auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "AUC (Area Under ROC Curve)",
  title = TeX("AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/ROC-AUC sensitivity Topology and Link Mutations (Bliss - Cascade 2.0)-2.png" width="80%" style="display: block; margin: auto;" />

### PR-AUC sensitivity (Bliss) {-}

Combine the $150$ simulation results (calibrated + proliferative)

```r
# Ensemble-wise
betas = seq(from = -10, to = 10, by = 0.1)

auc_values_ew = sapply(betas, function(beta) {
  pred_topolink_ew_bliss = pred_topolink_ew_bliss %>% mutate(combined_score = ss_score_150sim + beta * rand_score_150sim)
  res = pr.curve(scores.class0 = pred_topolink_ew_bliss %>% pull(combined_score) %>% (function(x) {-x}), 
    weights.class0 = pred_topolink_ew_bliss %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_ew = as_tibble(cbind(betas, auc_values_ew))

ggline(data = df_ew, x = "betas", y = "auc_values_ew", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("$\\beta$"), ylab = "AUC (Area Under PR Curve)",
  title = TeX("AUC sensitivity to $\\beta$ parameter: $calibrated + \\beta \\times random$"),
  color = my_palette[2], label = "betas", repel = TRUE,
  label.select = list(criteria = "`y` == max(`y`)"),
  font.label = list(size = 14)) + geom_vline(xintercept = 0) + grids()
```

<img src="index_files/figure-html/PR-AUC sensitivity Topology and Link Mutations (Bliss - Cascade 2.0)-1.png" width="80%" style="display: block; margin: auto;" />

```r
# Model-wise
weights = seq(from = 0, to = 1, by = 0.05)

auc_values_mw = sapply(weights, function(w) {
  pred_topolink_mw_bliss = pred_topolink_mw_bliss %>% 
    mutate(weighted_prob = (1 - w) * pred_topolink_mw_bliss$synergy_prob_ss_150sim + w * pred_topolink_mw_bliss$synergy_prob_rand_150sim)
  res = pr.curve(scores.class0 = pred_topolink_mw_bliss %>% pull(weighted_prob), 
    weights.class0 = pred_topolink_mw_bliss %>% pull(observed))
  auc_value = res$auc.davis.goadrich
})

df_mw = as_tibble(cbind(weights, auc_values_mw))

ggline(data = df_mw, x = "weights", y = "auc_values_mw", numeric.x.axis = TRUE,
  plot_type = "l", xlab = TeX("weight $w$"), ylab = "AUC (Area Under PR Curve)",
  title = TeX("AUC sensitivity to weighted average score: $(1-w) \\times prob_{ss} + w \\times prob_{rand}$"),
  color = my_palette[3]) + grids()
```

<img src="index_files/figure-html/PR-AUC sensitivity Topology and Link Mutations (Bliss - Cascade 2.0)-2.png" width="80%" style="display: block; margin: auto;" />

### ROC and PRC for best beta (HSA & Bliss) {-}

Use the best betas from the PR-AUC sensitivity graphs (separately for Bliss and HSA)!


```r
best_beta_hsa = -1
best_best_bliss = -1.3
pred_topolink_ew_hsa = pred_topolink_ew_hsa %>% mutate(best_score = ss_score_150sim + best_beta_hsa * rand_score_150sim)
pred_topolink_ew_bliss = pred_topolink_ew_bliss %>% mutate(best_score = ss_score_150sim + best_best_bliss * rand_score_150sim)

roc_best_res_hsa = get_roc_stats(df = pred_topolink_ew_hsa, pred_col = "best_score", label_col = "observed")
roc_best_res_bliss = get_roc_stats(df = pred_topolink_ew_bliss, pred_col = "best_score", label_col = "observed")

pr_best_res_hsa = pr.curve(scores.class0 = pred_topolink_ew_hsa %>% pull(best_score) %>% (function(x) {-x}), 
    weights.class0 = pred_topolink_ew_hsa %>% pull(observed), curve = TRUE, rand.compute = TRUE)
pr_best_res_bliss = pr.curve(scores.class0 = pred_topolink_ew_bliss %>% pull(best_score) %>% (function(x) {-x}), 
    weights.class0 = pred_topolink_ew_bliss %>% pull(observed), curve = TRUE)

# Plot best ROCs
plot(x = roc_best_res_hsa$roc_stats$FPR, y = roc_best_res_hsa$roc_stats$TPR,
  type = 'l', lwd = 3, col = my_palette[1], main = TeX('ROC curve (Ensemble-wise), $calibrated + \\beta \\times proliferative$'),
  xlab = 'False Positive Rate (FPR)', ylab = 'True Positive Rate (TPR)')
lines(x = roc_best_res_bliss$roc_stats$FPR, y = roc_best_res_bliss$roc_stats$TPR,
  lwd = 3, col = my_palette[3])
legend('bottomright', title = TeX('AUC ($\\beta_{HSA}$ = -1, $\\beta_{Bliss}$ = -1.3)'), 
  col = c(my_palette[1], my_palette[3]), pch = 19,
  legend = c(paste(round(roc_best_res_hsa$AUC, digits = 2), 'HSA (150 sim)'),
    paste(round(roc_best_res_bliss$AUC, digits = 2), 'Bliss (150 sim)')))
grid(lwd = 0.5)
abline(a = 0, b = 1, col = 'lightgrey', lty = 'dotdash', lwd = 1.2)

# Plot best PRCs
plot(pr_best_res_hsa, main = TeX('PR curve (Ensemble-wise), $calibrated + \\beta \\times proliferative$'),
  auc.main = FALSE, color = my_palette[1], rand.plot = TRUE)
plot(pr_best_res_bliss, add = TRUE, color = my_palette[3])
legend('topright', title = TeX('AUC ($\\beta_{HSA}$ = -1, $\\beta_{Bliss}$ = -1.3)'), col = c(my_palette[1], my_palette[3]), pch = 19,
  legend = c(paste(round(pr_best_res_hsa$auc.davis.goadrich, digits = 2), 'HSA (150 sim)'),
    paste(round(pr_best_res_bliss$auc.davis.goadrich, digits = 2), 'Bliss (150 sim)')))
grid(lwd = 0.5)
```

<img src="index_files/figure-html/ROC and PRC for best beta - Topology and Link Mutations (HSA and Bliss - Cascade 2.0)-1.png" width="50%" /><img src="index_files/figure-html/ROC and PRC for best beta - Topology and Link Mutations (HSA and Bliss - Cascade 2.0)-2.png" width="50%" />

## Fitness vs performance {-}

The idea here is to generate many training data files from the steady state as used in the simulations above for the AGS, where some of the nodes will have their states *flipped* to the opposite state ($0$ to $1$ and vice versa).
That way, we can train models to different steady states, ranging from ones that differ to just a few nodes states up to a steady state that is the complete *reversed* version of the one used in the simulations above.

Using the `gen_training_data.R` script, we first choose a few number of flips ($11$ flips) ranging from $1$ to $24$ (all nodes) in the steady state.
Then, for each such *flipping-nodes* value, we generated $20$ new steady states with a randomly chosen set of nodes whose value is going to flip.
Thus, in total, $205$ training data files were produced ($205 = 9 \times 20 + 24 + 1$, where from the $11$ number of flips, the one flip happens for every node and flipping all the nodes simultaneously happens once).

Running the script `run_druglogics_synergy_training.sh` from the `druglogics-synergy` repository root, we get the simulation results for each of these training data files.
The only difference in the cascade 2.0 configuration file was the number of simulations ($15$) for each training data file and the attractor tool used (`biolqm_stable_states`).

We now load the data from these simulations:

# Reproduce simulation results {-}

## ROC curves {-}

- Install `druglogics-synergy`
- Run the script `run_druglogics_synergy.sh` in the above repo using the configuration settings: 
  - `simulations: 50`
  - `attractor_tool: biolqm_stable_states`
  - `synergy_method: hsa` (also rerun the script chaning the synergy method `bliss`)

Thus you will get a directory per simulation and inside will be several result files.
To generate the ROC curves, we use the ensemble-wise and model-wise synergies found in each respective simulation.

## Random model results {-}

The CASCADE 1.0 and 2.0 `.sif` network files can be found at the directories `ags_cascade_1.0` and `ags_cascade_2.0` on the
`druglogics-synergy` repository.

Run the `abmlog` for the CASCADE 2.0 topology:
```
java -cp target/abmlog-1.5.0-jar-with-dependencies.jar eu.druglogics.abmlog.RandomBooleanModelGenerator --file=test/cascade_2_0.sif --num=3000
```

Next, prune the resulting models to only the ones that have 1 stable state ($1292$) using the simple bash script `process_models.sh` inside the generated `models` directory from `abmlog`.

```
cd pathTo/druglogics-synergy/ags_cascade_2.0
```

- Move the `models` dir inside the `ags_cascade_2.0` dir
- Use attractor_tool: `biolqm_stable_states` in the config file
- Use `synergy_method: hsa` or `synergy_method: bliss` (or run twice)
- Run drabme via `druglogics-synergy`:

```
java -cp ../target/synergy-1.2.0-jar-with-dependencies.jar eu.druglogics.drabme.Launcher --project=cascade_2.0_random_hsa --modelsDir=models --drugs=drugpanel --perturbations=perturbations --config=config --modeloutputs=modeloutputs
java -cp ../target/synergy-1.2.0-jar-with-dependencies.jar eu.druglogics.drabme.Launcher --project=cascade_2.0_random_bliss --modelsDir=models --drugs=drugpanel --perturbations=perturbations --config=config --modeloutputs=modeloutputs
```

The above procedure is the same for CASCADE 1.0. Changes:

- Network file is now the `cascade_1_0.sif`
- The `models` directory should be put inside the `ags_cascade_1.0` of `druglogics-synergy`
- The drabme command should be run with `--project=cascade_1.0_random_hsa` and `--project=cascade_1.0_random_bliss` respectively

# R session info {-}


```r
xfun::session_info()
```

```
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.4 LTS

Locale:
  LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
  LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
  LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
  LC_PAPER=en_US.UTF-8       LC_NAME=C                 
  LC_ADDRESS=C               LC_TELEPHONE=C            
  LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

Package version:
  assertthat_0.2.1        backports_1.1.6         base64enc_0.1.3        
  BH_1.72.0.3             bibtex_0.4.2.2          bookdown_0.18          
  broom_0.5.6             callr_3.4.3             Ckmeans.1d.dp_4.3.2    
  cli_2.0.2               clipr_0.7.0             codetools_0.2-16       
  colorspace_1.4-1        compiler_3.6.3          corrplot_0.84          
  cowplot_1.0.0           crayon_1.3.4            crosstalk_1.1.0.1      
  desc_1.2.0              digest_0.6.25           dplyr_0.8.5            
  DT_0.13                 ellipsis_0.3.0          emba_0.1.4             
  equatiomatic_0.0.0.9000 evaluate_0.14           fansi_0.4.1            
  farver_2.0.3            foreach_1.5.0           gbRd_0.4-11            
  generics_0.0.2          ggplot2_3.3.0           ggpubr_0.2.5           
  ggrepel_0.8.2           ggsci_2.9               ggsignif_0.6.0         
  glmnet_3.0-2            glue_1.4.0              graphics_3.6.3         
  grDevices_3.6.3         grid_3.6.3              gridExtra_2.3          
  gtable_0.3.0            highr_0.8               hms_0.5.3              
  htmltools_0.4.0         htmlwidgets_1.5.1       igraph_1.2.5           
  isoband_0.2.1           iterators_1.0.12        jsonlite_1.6.1         
  knitr_1.28              labeling_0.3            later_1.0.0            
  latex2exp_0.4.0         lattice_0.20-41         lazyeval_0.2.2         
  lifecycle_0.2.0         magrittr_1.5            MAMSE_0.2-1            
  markdown_1.1            MASS_7.3.51.5           Matrix_1.2-18          
  methods_3.6.3           mgcv_1.8.31             mime_0.9               
  munsell_0.5.0           nlme_3.1-145            pillar_1.4.3           
  pkgbuild_1.0.6          pkgconfig_2.0.3         pkgload_1.0.2          
  plogr_0.2.0             plyr_1.8.6              polynom_1.4.0          
  praise_1.0.0            prettyunits_1.1.1       processx_3.4.2         
  promises_1.1.0          PRROC_1.3.1             ps_1.3.2               
  purrr_0.3.3             R6_2.4.1                RColorBrewer_1.1-2     
  Rcpp_1.0.4.6            Rdpack_0.11-1           readr_1.3.1            
  reshape2_1.4.4          rje_1.10.15             rlang_0.4.5            
  rmarkdown_2.1           rprojroot_1.3.2         rstudioapi_0.11        
  scales_1.1.0            shape_1.4.4             splines_3.6.3          
  stats_3.6.3             stringi_1.4.6           stringr_1.4.0          
  testthat_2.3.2          tibble_3.0.0            tidyr_1.0.2            
  tidyselect_1.0.0        tinytex_0.21            tools_3.6.3            
  usefun_0.4.5            utf8_1.1.4              utils_3.6.3            
  vctrs_0.2.4             viridisLite_0.3.0       visNetwork_2.0.9       
  withr_2.1.2             xfun_0.12               yaml_2.2.1             
```

# References {-}
