#!/bin/bash

  for p in 4 6 8
  do
    for g in `seq 2500 2500 25000`
    do
    
      echo $p " " $g

      if [ $p == 4 ]
        then
          Rscript --vanilla alloSNPdrift.R 2 2 $g 100 20 \
          simAlloSNPdrift/alloSNP-p$p-c20-i100-g$g \
          simAlloSNP/alloSNP-p$p-c20-i100-freqs1.txt
      fi

      if [ $p == 6 ]
        then
          Rscript --vanilla alloSNPdrift.R 2 4 $g 100 20 \
          simAlloSNPdrift/alloSNP-p$p-c20-i100-g$g \
          simAlloSNP/alloSNP-p$p-c20-i100-freqs1.txt
      fi

      if [ $p == 8 ]
        then
          Rscript --vanilla alloSNPdrift.R 4 4 $g 100 20 \
          simAlloSNPdrift/alloSNP-p$p-c20-i100-g$g \
          simAlloSNP/alloSNP-p$p-c20-i100-freqs1.txt
      fi

    done
  done
