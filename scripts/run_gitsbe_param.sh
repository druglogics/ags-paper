#!/bin/bash

# Execute this script inside the `ags_cascade_2.0` directory of the
# `druglogics-synergy` repository (v1.2.1)

# config `simulations` was changed to 1500 to generate 4500 models
sed -i 's/simulations:.*/simulations:\t1500/' config

# link-operator mutations only: `1.2.1` version has already that configuration
sed -i 's/balance_mutations:.*/balance_mutations:\t3/' config
sed -i 's/topology_mutations:.*/topology_mutations:\t0/' config

echo "Running Gitsbe for 1500 simulations, fitting to steady state, link-operator mutations only"
start=`date +%s`
java -cp ../target/synergy-1.2.1-jar-with-dependencies.jar eu.druglogics.gitsbe.Launcher --project=gitsbe_link_only_cascade_2.0_ss --network=network.sif --trainingdata=training --config=config --drugs=drugpanel --modeloutputs=modeloutputs > /dev/null 2>&1
runtime=$(($(date +%s)-$start))
echo -e Execution Time: "$(($runtime / 60)) minutes\n"

# topology mutations only: no balance mutations, add topology ones
sed -i 's/balance_mutations:.*/balance_mutations:\t0/' config
sed -i 's/topology_mutations:.*/topology_mutations:\t10/' config

echo "Running Gitsbe for 1500 simulations, fitting to steady state, topology mutations only"
start=`date +%s`
java -cp ../target/synergy-1.2.1-jar-with-dependencies.jar eu.druglogics.gitsbe.Launcher --project=gitsbe_topology_only_cascade_2.0_ss --network=network.sif --trainingdata=training --config=config --drugs=drugpanel --modeloutputs=modeloutputs > /dev/null 2>&1
runtime=$(($(date +%s)-$start))
echo -e Execution Time: "$(($runtime / 60)) minutes\n"

# link-operator and topology mutations
sed -i 's/balance_mutations:.*/balance_mutations:\t3/' config
sed -i 's/topology_mutations:.*/topology_mutations:\t10/' config

echo "Running Gitsbe for 1500 simulations, fitting to steady state, both kind of mutations"
start=`date +%s`
java -cp ../target/synergy-1.2.1-jar-with-dependencies.jar eu.druglogics.gitsbe.Launcher --project=gitsbe_topo_and_link_cascade_2.0_ss --network=network.sif --trainingdata=training --config=config --drugs=drugpanel --modeloutputs=modeloutputs > /dev/null 2>&1
runtime=$(($(date +%s)-$start))
echo -e Execution Time: "$(($runtime / 60)) minutes\n"

