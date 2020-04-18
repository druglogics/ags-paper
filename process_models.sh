#!/bin/bash

# remove gitsbe files that have more or less than 1 stablestate

files=`ls | grep gitsbe`
for file in ${files}; do
  ss_num=`cat $file | grep stablestate | wc -l`
  if [ ${ss_num} != 1 ] 
  then
    rm $file
  fi
  # change modelname inside file to mimic gitsbe generated modelnames,
  # otherwise drabme will fail to recognise it!
  sed -i 's/modelname: cascade_2_0_/modelname: cascade_2_0_run_/' $file
done
