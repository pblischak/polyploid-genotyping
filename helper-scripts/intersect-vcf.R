#!/usr/bin/Rscript

start <- proc.time()

# Check to see if argparse, readr, and dplyr packages are installed.
# If not, install it.
if(!(c("argparse") %in% installed.packages()[,"Package"])){
  install.packages("argparse")
}

if(!(c("readr") %in% installed.packages()[,"Package"])){
  install.packages("readr")
}

if(!(c("dplyr") %in% installed.packages()[,"Package"])){
  install.packages("dplyr")
}

# import argparse library
suppressMessages(library(argparse))
library(readr)

# Set up cmd line arguments
parser <- ArgumentParser(description="Intersect two VCF files for shared variants.")
parser$add_argument("--vcf1", action="store", help="1st VCF file.")
parser$add_argument("--vcf2", action="store", help="2nd VCF file.")
parser$add_argument("--prefix", action="store", help="Prefix for output files.", default="shared")

# Assign parsed arguments to local variables
args     <- parser$parse_args()
vcf1     <- args$vcf1
vcf2     <- args$vcf2
prefix   <- args$prefix

# Read in each file as a data.frame
cat("\nReading in VCF files...\n")
vcf1_df <- read.table(vcf1, header=F)
vcf2_df <- read.table(vcf2, header=F)

# Find common variants using dplyr::inner_join
cat("Finding shared variants...\n")
var_df <- suppressMessages(dplyr::inner_join(vcf1_df[,c(1,2)],vcf2_df[,c(1,2)]))
cat("  ", nrow(var_df), " variants in common.\n", sep='')
write_tsv(var_df, paste(prefix,"-variants.txt", sep=''), col_names = F)

# Select the variants that are in var_df within each file and write
# to a new file with the PREFIX added.
cat("Writing shared variants to new VCF files...\n")
write_tsv(suppressMessages(dplyr::semi_join(vcf1_df, var_df)), paste(prefix, "-vcf1.vcf", sep=''), col_names = F)
write_tsv(suppressMessages(dplyr::semi_join(vcf2_df, var_df)), paste(prefix, "-vcf2.vcf", sep=''), col_names = F)

cat("Done.\n\n")

end <- proc.time() - start
cat("Elapsed time: ", end[3], " seconds.\n\n", sep='')
