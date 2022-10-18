#!/bin/bash
pref="../treewidth_solver_input/"
for filename in ../treewidth_solver_input/*.txt; do
      # echo ${filename}
      echo "./flow_cutter_pace17 < ${filename} > ../treewidth_solver_output/${filename/#$pref}"
      ./flow_cutter_pace17 < ${filename} > ../treewidth_solver_output/${filename/#$pref} &
      sleep 30
      kill -SIGINT "$!"
done
