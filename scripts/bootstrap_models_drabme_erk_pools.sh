#!/bin/bash

# Execute this script inside the `ags_cascade_2.0` directory of the 
# `druglogics-synergy` repository (v1.2.0)

# Use this script to bootstrap `batch_size` amount of gitsbe models from
# a specific `models_dir` directory and provide them as input to `drabme`.
# This process will be repeated `batches` times.

batches=35
batch_size=300

# change drabme config option (use bliss synergy assessement)                   
# all other options to default values in v1.2.0                                 
sed -i 's/synergy_method:.*/synergy_method:\tbliss/' config

# from the Zenodo dataset http://tiny.cc/ags-paper-zenodo, file `erk_perf_investigation.tar.gz`,
# copy the `erk_perf_investigation/erk_active_pool` and the
# `erk_perf_investigation/erk_inhibited_pool` directories,
# into the `ags_cascade_2.0` directory of the `druglogics-synergy` repository

echo Simulations using files with ERK ACTIVE
models_dir="erk_active_pool"

for batch in $( seq 1 $batches )
do
  echo "Batch No.${batch}"
  sample_models_dir=${models_dir}_${batch}

  # create directory to store the sampled models
  mkdir $sample_models_dir

  # copy sampled models
  ls ${models_dir} | sort -R | tail -$batch_size | while read file; do
    cp ${models_dir}/$file $sample_models_dir
  done

  # run Drabme
  start=`date +%s`
  java -cp ../target/synergy-1.2.0-jar-with-dependencies.jar eu.druglogics.drabme.Launcher --project=cascade_2.0_${models_dir}_bliss_batch_${batch} --modelsDir=${sample_models_dir} --drugs=drugpanel --perturbations=perturbations --config=config --modeloutputs=modeloutputs > /dev/null 2>&1
  runtime=$(($(date +%s)-$start))
  echo -e Execution Time: "$(($runtime / 60)) minutes\n"
done

echo Simulations using files with ERK INHIBITED
models_dir="erk_inhibited_pool"

for batch in $( seq 1 $batches )                                                
do                                                                              
  echo "Batch No.${batch}"                                                      
  sample_models_dir=${models_dir}_${batch}                                      
                                                                                
  # create directory to store the sampled models                                
  mkdir $sample_models_dir                                                      
                                                                                
  # copy sampled models                                                         
  ls ${models_dir} | sort -R | tail -$batch_size | while read file; do          
    cp ${models_dir}/$file $sample_models_dir                                   
  done                                                                          
                                                                                
  # run Drabme                                                                  
  start=`date +%s`                                                              
  java -cp ../target/synergy-1.2.0-jar-with-dependencies.jar eu.druglogics.drabme.Launcher --project=cascade_2.0_${models_dir}_bliss_batch_${batch} --modelsDir=${sample_models_dir} --drugs=drugpanel --perturbations=perturbations --config=config --modeloutputs=modeloutputs > /dev/null 2>&1
  runtime=$(($(date +%s)-$start))                                               
  echo -e Execution Time: "$(($runtime / 60)) minutes\n"                        
done

