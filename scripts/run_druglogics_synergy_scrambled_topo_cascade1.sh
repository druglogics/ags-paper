#!/bin/bash

# Run this script from the `druglogics-synergy` repository root (v1.2.0)
# The `scrambled_topologies_cascade1` directory is produced via the
# `gen_scrambled_topologies_cascade1.R` script and should also be placed in
# the repository root

# The result files/directories are stored in Zenodo [TOADD],
# file `synergy-res-scrambled-topo-cascade1.tar.gz`

training_files=`ls scrambled_topologies_cascade1`
files_num=`ls scrambled_topologies_cascade1 | wc -l`
count=0

# change drabme config option (use bliss synergy assessement)
# all other options to default values in v1.2.0
sed -i 's/synergy_method:.*/synergy_method:\tbliss/' ags_cascade_1.0/config

for file in ${training_files}; do
  count=$((count + 1))
  echo Using topology file No. $count/$files_num: $file
  cat scrambled_topologies_cascade1/$file > ags_cascade_1.0/network.sif

  # strip .sif extension for project naming
  file_name=$(echo "$file" | sed 's/\.[^.]*$//')

  start=`date +%s`

  # run calibrated models simulations
  cat ags_cascade_1.0/steadystate > ags_cascade_1.0/training
  java -cp target/synergy-1.2.0-jar-with-dependencies.jar eu.druglogics.synergy.Launcher --inputDir=ags_cascade_1.0 --project=${file_name}_ss > /dev/null 2>&1

  # run random proliferative models simulations
  cat ags_cascade_1.0/random_train > ags_cascade_1.0/training
  java -cp target/synergy-1.2.0-jar-with-dependencies.jar eu.druglogics.synergy.Launcher --inputDir=ags_cascade_1.0 --project=${file_name}_random > /dev/null 2>&1

  runtime=$(($(date +%s)-$start))
  echo Execution Time: "$(($runtime / 60)) minutes"
done

