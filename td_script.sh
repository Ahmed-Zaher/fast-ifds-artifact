#!/bin/bash
pref="../treedepth_solver_input/"
for filename in ../treedepth_solver_input/*.txt; do
      echo ${filename}
      ./flow_cutter_parallel_pace20 < ${filename} > ../treedepth_solver_output/${filename/#$pref} &
      sleep 30
      kill -SIGINT "$!"
done
