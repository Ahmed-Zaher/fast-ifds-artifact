#!/bin/bash
pref="../treewidth_solver_input/"
for filename in ../treewidth_solver_input/*.txt; do
      # echo ${filename}
      echo "./tw-heuristic < ${filename} > ../treewidth_solver_output/${filename/#$pref}"
      ./tw-heuristic < ${filename} > ../treewidth_solver_output/${filename/#$pref} &
      sleep 30
      kill -SIGINT "$!"
done
