#!/bin/bash

# Execute this script inside the `ags_cascade_2.0` directory of the
# `druglogics-synergy` repository (v1.2.1)

# Use this script to bootstrap `batch_size` amount of gitsbe models from
# a specific `models_dir` directory and provide them as input to `drabme`.
# This process will be repeated `batches` times.

batches=20
batch_size=100

# the `models` directory can for example be created with the `run_gitsbe_random.sh` script
models_dir="models"

for batch in $( seq 1 $batches )
do
  echo "Batch No.${batch}"
  sample_models_dir=models_batch_${batch}

  # create directory to store the sampled models
  mkdir $sample_models_dir

  # copy sampled models
  ls ${models_dir} | sort -R | tail -$batch_size | while read file; do
    cp ${models_dir}/$file $sample_models_dir
  done

  # run Drabme
  start=`date +%s`
  java -cp ../target/synergy-1.2.1-jar-with-dependencies.jar eu.druglogics.drabme.Launcher --project=cascade_2.0_rand_prolif_bliss_batch_${batch} --modelsDir=${sample_models_dir} --drugs=drugpanel --perturbations=perturbations --config=config --modeloutputs=modeloutputs > /dev/null 2>&1
  runtime=$(($(date +%s)-$start))
  echo -e Execution Time: "$(($runtime / 60)) minutes\n"
done

