#!/bin/bash

# Run this script from the `druglogics-synergy` repository root (v1.2.1)
# The `training-data-files` are produced via the `gen_training_data.R` script
# See Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `training-data-files.tar.gz`
training_files=`ls training-data-files`
files_num=`ls training-data-files | wc -l`
count=0

# change `ags_cascade_1.0` to `ags_cascade_2.0` depending on which topology you want to use

# change number of simulations
sed -i 's/simulations:.*/simulations:\t50/' ags_cascade_1.0/config # default value is 50
# sed -i 's/simulations:.*/simulations:\t20/' ags_cascade_2.0/config # 20 for balance mutation, 50 for topology mutations

# change drabme config option (use bliss synergy assessement)
sed -i 's/synergy_method:.*/synergy_method:\tbliss/' ags_cascade_1.0/config
# sed -i 's/synergy_method:.*/synergy_method:\tbliss/' ags_cascade_2.0/config

for file in ${training_files}; do
  count=$((count + 1))
  echo Using flipped training file No. $count/$files_num: $file
  cat training-data-files/$file > ags_cascade_1.0/training
  #cat training-data-files/$file > ags_cascade_2.0/training

  start=`date +%s`
  java -cp target/synergy-1.2.1-jar-with-dependencies.jar eu.druglogics.synergy.Launcher --inputDir=ags_cascade_1.0 --project=$file > /dev/null 2>&1 # for CASCADE 2.0: `--inputDir=ags_cascade_2.0`
  runtime=$(($(date +%s)-$start))
  echo Execution Time: "$(($runtime / 60)) minutes"
done
