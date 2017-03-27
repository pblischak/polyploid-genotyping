## `Rcode/`

These are the R, C++, and bash scripts that were used to conduct the simulation study for our manuscript.
To run them, you will need to have the R package `Rcpp` installed. They take a pretty long time to run (several hours),
so be aware when trying to redo them.

### `diseq` simulations

Execute the bash script `simDiseq.sh`:

```bash
# run simDiseq bash script
bash simDiseq.sh

# change into the simDiseq/ directory and analyze the data with ebg
cd simDiseq/
bash runDiseq.sh      # Analyze with the diseq model
bash runDiseqGATK.sh  # Analyze with the gatk model
bash runDiseqHWE.sh   # Analyze with the hwe model
```

### `alloSNP` simulations

Execute the bash script `simAlloSNP.sh`:

```bash
# run simAlloSNP bash script
bash simAlloSNP.sh

# change into the simAlloSNP/ directory and analyze the data with ebg
cd simAlloSNP/
bash runAlloSNP.sh       # Analyze with the alloSNP model
bash runAlloSNPbrent.sh  # Analyze with the alloSNP model EM+Brent algorithm
bash runAlloSNPhwe.sh    # Analyze with the hwe model
bash runAlloSNPgatk.sh   # Analyze with the gatk model
```
