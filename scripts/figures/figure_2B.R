library(readr)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggpubr)

# define a function to read a *summary.txt output gitsbe file
read_summary_file = function(file_name) {
  lines = readr::read_lines(file = file_name, skip = 5, skip_empty_rows = TRUE)

  data_list = list()
  index = 1

  gen_fit_list = list()
  gen_index = 1
  for (line_index in 1:length(lines)) {
    line = lines[line_index]
    if (stringr::str_detect(string = line, pattern = "Simulation")) {
      data_list[[index]] = dplyr::bind_cols(gen_fit_list)
      index = index + 1

      gen_fit_list = list()
      gen_index = 1
    } else { # read fitness values
      gen_fit_list[[gen_index]] = tibble::as_tibble_col(as.numeric(unlist(strsplit(line, split = '\t'))), column_name = paste0(gen_index))
      gen_index = gen_index + 1
    }
  }

  # add the last simulation's values
  data_list[[index]] = dplyr::bind_cols(gen_fit_list)

  return(data_list)
}

# input summary file
# generated using CASCADE 1.0 topology, with a configuration of 1000 simulations,
# `bootstrap_mutations_factor` equal to 1 and hsa as `synergy_method`
fitness_summary_file = "results/link-only/cascade_1.0_ss_1000sim_fixpoints_hsa_summary.txt"

# `fit_res` is a list of tibbles
# Each tibble has the fitness results of a simulation
# Rows represent the models and columns are the generations
fit_res = read_summary_file(file_name = fitness_summary_file)

# average fitness + standard deviation per generation across all 1000 simulations
# `avg_fit` is a tibble with rows the number of simulations and
# columns the generations. Each value in a (sim,gen) cell is the average
# fitness of models in that particular (sim,gen) combination
avg_fit = do.call(dplyr::bind_rows, sapply(fit_res, colMeans))
avg_fit_long = avg_fit %>%
  tidyr::pivot_longer(cols = everything()) %>%
  mutate(name = as.integer(name))

set1_cols = RColorBrewer::brewer.pal(n = 9, name = "Set1")

ggline(data = avg_fit_long, x = "name", y = "value", color = set1_cols[2],
  add = "mean_sd", add.params = list(color = "black"), ylim = c(0, 1),
  main = "Fitness Evolution across Generations",
  xlab = "Generations", ylab = "Fitness") +
  theme(plot.title = element_text(hjust = 0.5)) + grids()
ggsave(filename = 'scripts/figures/figure_2B.png', dpi = 400, width = 7, height = 5)
