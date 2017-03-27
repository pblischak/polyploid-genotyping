#!/bin/bash

for i in 25 50 100
do
  for p in 4 6 8
  do
    for c in 2 5 10 20 30 40
    do
          #Rscript --vanilla alloSNP.R 2 2 $i $c simAlloSNP/alloSNP-p$p-c$c-i$i
        ebg gatk -n $i -l 10000 -p $p \
            -t alloSNP-p$p-c$c-i$i-tot.txt \
            -r alloSNP-p$p-c$c-i$i-ref.txt \
            -e alloSNP-p$p-c$c-i$i-err.txt \
            --prefix gatk-results/p$p-i$i-c$c

    done
  done
done
