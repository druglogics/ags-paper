#!/bin/bash

# remove gitsbe files that have more or less than 1 stablestate

files=`ls | grep gitsbe`
for file in ${files}; do
  ss_num=`cat $file | grep stablestate | wc -l`
  if [ ${ss_num} != 1 ] 
  then
    rm $file
  else
    # change modelname inside file to mimic gitsbe generated modelnames,
    # otherwise drabme will fail to recognise it!
    # Remember to change 'cascade_1_0' to 'cascade_2_0' if appropriate
    sed -i 's/modelname: cascade_1_0_/modelname: cascade_1_0_run_/' $file
  fi
done
