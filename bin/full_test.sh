#!/bin/bash

server='http://localhost:8080'
for paradigm_in_degree in 5 7 ; do
  for ceres_threshold in -0.8 -0.6 -1 ; do
    for use_steady_state in true false ; do
      for use_observed_synergies in true false ; do
        for obs_method in avg_bliss GS; do
          for dir in AGS.slow Topology20171116; do
            for cell_line in AGS SW620 Colo205 SF295; do
                cmd="drbbt task SINTEF ROC --load_inputs /home/miguelvg/git/workflows/SINTEF/examples/ROC/$dir --ceres_threshold $ceres_threshold --use_steady_state $use_steady_state --use_observed_synergies $use_observed_synergies --obs_method $obs_method --cell_line $cell_line --jn $cell_line -pf --log 1"
                echo -n "RUN: " >>  log
                echo $cmd >> log
                file=`$cmd 2>> log`
                echo -n "RUN: "
                echo $cmd
                echo -n "FILE: "
                echo $file
                echo -n "FILE: " >> log
                echo $file >> log
                basename=`basename $file`
                url="$server/SINTEF/ROC/$basename"
                echo -n "URL: "
                echo $url
                echo -n "URL: " >> log
                echo $url >> log
            done
          done
        done
      done
    done
  done
done
