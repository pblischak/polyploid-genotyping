#!/bin/bash

# Run analysis on all simulated data sets for the diseq model

for i in 25 50 100
do
  for p in 4 6 8
  do
    for c in 2 5 10 20 30
    do
      for f in 10 25 50 75 90
      do

        #echo $(bc <<< "scale=2; $f / 100")
        #Rscript --vanilla diseq.R $(bc <<< "scale=2; $f / 100") $p $i $c simDiseq/diseq-F$f-p$p-i$i-c$c
        ebg hwe -n $i -l 10000 -p $p \
            -t diseq-F$f-p$p-i$i-c$c-tot.txt \
            -r diseq-F$f-p$p-i$i-c$c-ref.txt \
            -e diseq-F$f-p$p-i$i-c$c-err.txt \
            --iters 1000 \
            --prefix hwe-results/F$f-p$p-i$i-c$c

      done
    done
  done
done
