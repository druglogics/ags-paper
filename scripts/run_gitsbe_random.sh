#!/bin/bash

# Execute this script inside the `ags_cascade_2.0` directory of the 
# `druglogics-synergy` repository (v1.2.0)

# change to random training data and 1000 simulations
cat random_train > training
sed -i 's/simulations:.*/simulations:\t1000/' config

echo "Running Gitsbe for 1000 simulations, fitting to a random proliferation profile, models with link-operator mutations"
start=`date +%s`
java -cp ../target/synergy-1.2.0-jar-with-dependencies.jar eu.druglogics.gitsbe.Launcher --project=cascade_2.0_rand_1000sim --network=network.sif --trainingdata=training --config=config --drugs=drugpanel --modeloutputs=modeloutputs > /dev/null 2>&1
runtime=$(($(date +%s)-$start))
echo -e Execution Time: "$(($runtime / 60)) minutes\n"

