## `ebg/`

`ebg` is a C++ program for estimating genotypes from high throughput sequencing data and works on both diploids and polyploids.
Input data include (1) a matrix of total read counts mapping to each site for each individual,
(2) a matrix of reference read counts mapping to each site for each individual, and
(3) a vector of sequencing error values for each locus.
Examples of each of these files can be found in the `data/` folder in the main GitHub repo and a walkthrough for how to analyze them is in the associated GitHub pages website.
Missing data are also allowed and should be given a value of `-9`. These counts and error values are usually found in VCF files as the allele depth (AD) field for genotypes and the QUAL field for each site.

There is a `Makefile` that will compile the `ebg` executable. Simply type:

```bash
# Compile executable
make

# Will cp ebg to /usr/local/bin
sudo make install
```

`ebg` is run from the command line and accepts parameter values as flags. To see a help message that displays the available models, type:

```bash
ebg -h
```

For help with any of the models implemented in `ebg`, just type one of the following:

```bash
# Help for the Hardy Weinberg model
ebg hwe -h

# Help for the disequilibrium model
ebg diseq -h

# Help for the allopolyploid subgenome model
ebg alloSNP -h
```

### `hwe`

The `hwe` model estimates genotypes assuming Hardy Weinberg equilibrium. You can run the model as follows:

```bash
ebg hwe -p <ploidy> \
  -n <num-individuals> \
  -l <num-loci> \
  -t <total-reads> \
  -r <ref-reads> \
  -e <error-vals>
```

### `diseq`

The `diseq` model estimates genotypes assuming Hardy Weinberg equilibrium. You can run the model as follows:

```bash
ebg hwe -p <ploidy> \
  -n <num-individuals> \
  -l <num-loci> \
  -t <total-reads> \
  -r <ref-reads> \
  -e <error-vals>
```

### `alloSNP`

The `alloSNP` model estimates genotypes within the subgenomes of an allopolyploid using a reference panel of allele frequencies estimated separately for one of the parents. You can run the model as follows:

```bash
ebg hwe -f <ref-panel-file>
  -n <num-individuals> \
  -l <num-loci> \
  -p1 <ploidy1> \
  -p2 <ploidy2>
  -t <total-reads> \
  -r <ref-reads> \
  -e <error-vals>
```

An option that is unique to the `alloSNP` model is the `--brent` flag, which is a logical flag indicating that Brent's method should be used if the algorithm hasn't converged after the specified number of EM iterations. It is usually a good idea to use it whenever you run an analysis with this model.

### Additional settings

These options can also be set for any of the models to control aspects of the EM algorithm
and how the output is written to file.

 - `--prefix <outfile-prefix>` : use the specified string as a prefix for the output files.
 - `--tol <tolerance>` : this is the tolerance used for Brent's method.
 - `--iters <number-iterations>` : the number of EM iterations to run.
 - `--stop <stop-criterion>` : this convergence criterion to stop updating a parameter
 - `--quiet` : a logical flag that suppressing the printing of output to stdout.
