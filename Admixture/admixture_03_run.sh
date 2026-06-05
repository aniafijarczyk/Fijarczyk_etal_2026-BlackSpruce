#!/bin/bash

#DIR=/home/software/admixture_linux-1.3.0

# chromosome names must be human (numeric)
bed=admixture_01_input.bed

#for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
for K in 2 3
  do
  echo ${K}
  admixture --cv ${bed} $K | tee log${K}.out
  done



