Scripts for filtering a VCF file and converting to read count matrices using
[VCFtools](https://vcftools.github.io/index.html) and a simple perl script (`read-counts-from-vcf.pl`).
The perl script takes input from `stdin` and is run with the following options passed to it (in order):

```perl
perl read-counts-from-vcf.pl <total-reads-file> <alt-reads-file> <allele-depth-position> <min-depth-filter>
```

 - `total-reads-file`: the name to be given to the total read count file.
 - `alt-reads-file`: the name to be given to the alternative read count file.
 - `allele-depth-position`: the position of the allele depth information in the genotype field for the VCF file.
 - `min-depth-filter`: the minimum number of reads for including read count data.

A typical use case would be to filter a VCF file for biallelic sites only, as well as applying whatever depth and quality filters
you would like to add. Using the `--recode` and `--stdout` flags with VCFtools will generate a new VCF file that is
printed to `stdout` that we can then pipe (`|`) into the perl script.

```bash
vcftools --gzvcf input.vcf.gz --max-alleles 2 --min-alleles 2 -minDP 5 --max-missing 0.5 \
         --recode --stdout | perl read-counts-from-vcf.pl tot.txt alt.txt 2 5
```

This command takes a gzipped VCF file, filters for biallelic sites, includes individuals with less than 50% missing data,
applies a minimum depth criterion and then prints a new VCF file to `stdout`. This new VCF file is piped as input to the perl
script, which prints the total and alternative read count depths for all individuals and sites to the files `tot.txt` and `alt.txt`,
respectively. The `2` tells the script that the allele depth (AD) field is the second entry in the genotype information encoded within the VCF file.
The last number is another read depth filter for the AD field, which ensures that all included sites have a minimum total number of reads greater than or equal to `5`.
