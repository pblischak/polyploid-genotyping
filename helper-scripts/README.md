These are the Python, R, and Bash scripts that were used to process the data files that
were generated as part of our manuscript. They mostly involve filtering and processing VCF files.
The Bash scripts typically also call external programs, so make sure that you have them installed:

 - [VCFtools](https://vcftools.github.io/index.html)
 - [BWA](http://bio-bwa.sourceforge.net/)
 - [SAMtools](http://www.htslib.org/doc/samtools.html)
 - [GATK](https://software.broadinstitute.org/gatk/)
 - [Picard](https://broadinstitute.github.io/picard/)

We used Python v2.7.13 and R v3.3.2.

### Scripts

The Bash scripts are included here as examples to show the commands that we used to process our data. They may not work for you "out of the box". The Python and R scripts should work for processing VCF files and all use a command line parsing library (argparse). To see the script options, run it with a `-h` flag.

**Analyses for _Andropogon gerargii_**

 - `andropogon-read-counts.sh`: Bash scripts that wraps a call to VCFtools and the Python script `read-counts-from-vcf.pl`.
 - `read-counts-from-vcf.pl`: Old Perl script that extracts the allele depth (AD) field from a VCF file and writes a file for the total reads and the number of alternative allele reads.
 - `filter-inds.R`: R script to filter sites output by `read-counts-from-vcf.pl` based on the number of missing individuals.
 ```
 Rscript filter-inds.R <tot-reads> <alt-reads> <%-missing> <transpose?> <missing-string>
 ```

**Analyses for _Betula_ species using GATK**

 - `index.sh`: Bash script to index the reference genome of *Betula nana* using bwa, SAMtools, and Picard.
 - `map_bwa.sh`: Bash script that loops through fastq files to map to the *Betaula* reference genome using BWA, followed by processing with SAMtools.
 - `add_read_groups.sh`: Bash script to add read group information with Picard.
 - `run_genotyper.sh`: Bash script to run the GATK UnifiedGenotyper on BAM files with read group information added by Picard.
 - `filter.vcf.R`: R script to filter a VCF file based on variant quality (QUAL), biallelic SNPs, sequencing read depth per genotype, and amount of missing data per site. Replaces the use of VCFtools and the `filter-inds.R` script.
 - `intersect-vcf.R`: Finds the shared variants between two VCF files and prints them to new files.
 - `read-counts-from-vcf.py`: Python script that replaces the old Perl version with a better interface and prints the output files in the correct orientation (doesn't need to be transposed).
 - `gt-from-vcf.py`: Python script to extract the count of alternative alleles in the GT field for each individual at each site.
 - `run_mpileup.sh`: Bash script to generate a pileup file from input BAM files.
 - `per-locus-err.py`: Extract the average PHRED-scaled base qualities for each site from a pileup file (default=0.01).
