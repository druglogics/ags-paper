#!/bin/bash

# remove gitsbe files that have more or less than 1 stablestate

files=`ls | grep gitsbe`
for file in ${files}; do
  ss_num=`cat $file | grep stablestate | wc -l`
  if [ ${ss_num} != 1 ] 
  then
    rm $file
  else
  fi
done
