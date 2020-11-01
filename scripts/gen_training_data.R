########################################
# Generate AGS training data files for #
# Fitness vs Performance Investigation #
########################################
library(dplyr)
library(tibble)

# the generated training data files are in the Zenodo dataset [TOADD],
# file `training-data-files.tar.gz`

# specify the output directory (and make sure it exists!)
output_dir = '/home/john/tmp/ags-paper/training-data-files'

# Read the AGS steady state and reformat it to a data.frame
steady_state_file = 'data/steadystate'
lines = readLines(steady_state_file)
steady_state = unlist(strsplit(x = lines[8], split = "\t"))

res_list = list()
for(el in steady_state) {
  res = unlist(strsplit(el, split = ":"))
  res_list[[el]] = c(res[1], res[2])
}

df = as.data.frame(do.call(rbind, res_list), stringsAsFactors = FALSE)
df = remove_rownames(df)
colnames(df) = c("nodes", "state")
df$state = as.integer(df$state)
df = tibble::as_tibble(df)

# Generate flipped steady states

set.seed(2)
num_of_flips = c(1,3,6,9,12,15,17,19,21,22,24)
random_states_num = 20
gold_state = df$state

flipped_states = list()
for(flip_num in num_of_flips) {
  print(flip_num)

  if (flip_num == 1) { # flip once every node
    for (flip_index in 1:nrow(df)) {
      flipped_state = gold_state
      flipped_state[flip_index] = !flipped_state[flip_index]
      flipped_states[[paste0(flip_num,"_flips_",flip_index)]] = flipped_state
    }
    next
  }

  if (flip_num == nrow(df)) { # flip all nodes at once!
    flipped_state = gold_state
    flipped_state = as.integer(!flipped_state)
    flipped_states[[paste0(flip_num,"_flips_1")]] = flipped_state
    next
  }

  for(state_index in 1:random_states_num) {
    flip_index_vec = sample(x = 1:nrow(df), size = flip_num)
    flipped_state = gold_state
    flipped_state[flip_index_vec] = !flipped_state[flip_index_vec]
    flipped_states[[paste0(flip_num,"_flips_",state_index)]] = flipped_state
  }
}

# Write flipped steady states to files
flipped_state_lines = lines

for (index in 1:length(flipped_states)) {
  flipped_state = flipped_states[[index]]
  state_name    = names(flipped_states)[index]

  flipped_state_lines[8] = paste0(df$nodes, ":", flipped_state, collapse = "\t")
  file_name = paste0(output_dir, "/training_", state_name)
  writeLines(text = flipped_state_lines, con = file_name)
}

