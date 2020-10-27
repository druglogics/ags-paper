library(dplyr)
library(tibble)
library(readr)
library(ggplot2)
library(usefun)

# the generated topology files are in the Zenodo dataset [TOADD],
# file `scrambled_topologies_cascade_1.0.tar.gz`

# specify the output directory
output_dir = '/home/john/tmp/ags-paper/scrambled_topologies_cascade1/'

# CASCADE 1.0
edge_tbl = readr::read_delim(file = 'https://raw.githubusercontent.com/druglogics/cascade/master/cascade_1.0.sif', delim = "\t", col_names = c('source', 'effect', 'target'), col_types = "ccc")

sources = edge_tbl %>% pull(source)
targets = edge_tbl %>% pull(target)
effects = edge_tbl %>% pull(effect)
edge_annot = edge_tbl %>%
  mutate(edge = paste(source, effect, target)) %>% pull(edge)

# for reproducibility
set.seed(42)

# expected similarity scores between each generated topology and CASCADE 1.0
exp_sims = c(seq(from = 0, to = 0.95, by = 0.05), 0.96, 0.98)

# total generated topologies: length(exp_sims) * 10 = 220 (x4)
count = 1
for (exp_sim in exp_sims) {
  for(i in 1:10) {
    # sources
    permut_sources = usefun::partial_permut(x = sources, exp_sim)
    # edge_annot_src_permut = paste(permut_sources, effects, targets)
    # src_sim = sum(edge_annot_src_permut %in% edge_annot)/length(edge_annot)

    # effects
    permut_effects = effects
    eff_indexes = sample(1:length(effects), size = round((1-exp_sim)*length(effects)))

    # change activation to inhibitions and vice-versa
    for (index in eff_indexes) {
      permut_effects[index] = ifelse(effects[index] == '->', '-|', '->')
    }
    # edge_annot_eff_permut = paste(sources, permut_effects, targets)
    # eff_sim = sum(edge_annot_eff_permut %in% edge_annot)/length(edge_annot)

    # targets
    permut_targets = usefun::partial_permut(x = targets, exp_sim)
    # edge_annot_trg_permut = paste(sources, effects, permut_targets)
    # trg_sim = sum(edge_annot_trg_permut %in% edge_annot)/length(edge_annot)

    # permute all (sources, effects and targets!)
    # edge_annot_all_permut = paste(permut_sources, permut_effects, permut_targets)
    # all_sim = sum(edge_annot_all_permut %in% edge_annot)/length(edge_annot)

    # Print similarity scores
    # print(paste0(exp_sim, " : ", src_sim, ",", trg_sim, ",", eff_sim))
    # print(paste0(exp_sim, " : ", all_sim))

    # save .sif files
    # sources scrambled
    readr::write_tsv(x = edge_tbl %>% mutate(source = permut_sources),
      file = paste0(output_dir, "/cascade1_src_permut_", count, ".sif"),
      col_names = FALSE)
    # targets scrambled
    readr::write_tsv(x = edge_tbl %>% mutate(target = permut_targets),
      file = paste0(output_dir, "/cascade1_trg_permut_", count, ".sif"),
      col_names = FALSE)
    # effects scrambled
    readr::write_tsv(x = edge_tbl %>% mutate(effect = permut_effects),
      file = paste0(output_dir, "/cascade1_eff_permut_", count, ".sif"),
      col_names = FALSE)
    # all 3 columns scrambled!
    readr::write_tsv(x = edge_tbl %>% mutate(source = permut_sources,
      effect = permut_effects, target = permut_targets),
      file = paste0(output_dir, "/cascade1_all_permut_", count, ".sif"),
      col_names = FALSE)
    count = count + 1
  }
}
