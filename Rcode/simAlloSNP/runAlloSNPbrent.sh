#!/bin/bash

for i in 25 50 100
do
  for p in 4 6 8
  do
    for c in 2 5 10 20 30 40
    do

      if [ $p == 4 ]
        then
          #Rscript --vanilla alloSNP.R 2 2 $i $c simAlloSNP/alloSNP-p$p-c$c-i$i
          ebg alloSNP -n $i -l 10000 -p1 2 -p2 2 \
          	  -f alloSNP-p$p-c$c-i$i-freqs1.txt \
              -t alloSNP-p$p-c$c-i$i-tot.txt\
              -r alloSNP-p$p-c$c-i$i-ref.txt\
              -e alloSNP-p$p-c$c-i$i-err.txt\
              --iters 100 --brent \
              --prefix brent-results/p$p-i$i-c$c
      fi

      if [ $p == 6 ]
        then
          #Rscript --vanilla alloSNP.R 2 4 $i $c simAlloSNP/alloSNP-p$p-c$c-i$i
          ebg alloSNP -n $i -l 10000 -p1 2 -p2 4 \
              -f alloSNP-p$p-c$c-i$i-freqs1.txt \
              -t alloSNP-p$p-c$c-i$i-tot.txt\
              -r alloSNP-p$p-c$c-i$i-ref.txt\
              -e alloSNP-p$p-c$c-i$i-err.txt\
              --iters 100 --brent \
              --prefix brent-results/p$p-i$i-c$c
      fi

      if [ $p == 8 ]
        then
          #Rscript --vanilla alloSNP.R 4 4 $i $c simAlloSNP/alloSNP-p$p-c$c-i$i
          ebg alloSNP -n $i -l 10000 -p1 4 -p2 4 \
              -f alloSNP-p$p-c$c-i$i-freqs1.txt \
              -t alloSNP-p$p-c$c-i$i-tot.txt\
              -r alloSNP-p$p-c$c-i$i-ref.txt\
              -e alloSNP-p$p-c$c-i$i-err.txt\
              --iters 100 --brent \
              --prefix brent-results/p$p-i$i-c$c
      fi

    done
  done
done
