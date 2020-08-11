#!/bin/bash

training_files=`ls ags_cascade_2.0/training-data-files` 
files_num=`ls ags_cascade_2.0/training-data-files | wc -l`
count=0

for file in ${training_files}; do
  count=$((count + 1))
  echo Using flipped training file No. $count/$files_num: $file
  cat ags_cascade_2.0/training-data-files/$file > ags_cascade_2.0/training
  
  start=`date +%s`
  java -cp target/synergy-1.2.0-jar-with-dependencies.jar eu.druglogics.synergy.Launcher --inputDir=ags_cascade_2.0 --project=$file > /dev/null 2>&1
  runtime=$(($(date +%s)-$start))
  echo Execution Time: "$(($runtime / 60)) minutes"
done
