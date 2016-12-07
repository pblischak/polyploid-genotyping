#!/bin/bash

# Simulate all data sets for the diseq model

for i in 25 50 100
do
  for p in 6 8
  do
    for c in 5 10 20
    do
      for f in 10 25 50 75 90
      do

        #echo $(bc <<< "scale=2; $f / 100")
        Rscript --vanilla diseq.R $(bc <<< "scale=2; $f / 100") $p $i $c simDiseq/diseq-F$f-p$p-i$i-c$c

      done
    done
  done
done
