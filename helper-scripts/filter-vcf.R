#!/usr/bin/Rscript

start <- proc.time()

# Check to see if argparse, readr packages are installed.
# If not, install it.
if(!(c("argparse") %in% installed.packages()[,"Package"])){
  install.packages("argparse")
}

if(!(c("readr") %in% installed.packages()[,"Package"])){
  install.packages("readr")
}

# define function to ceck missing data
check_missing <- function(x){
  res <- rep(NA,length(x))
  for(r in 1:length(res)){
    res[r] <- substr(as.character(x[r]),1,1)=="."
  }
  return(res)
}

# define function to extract DP field from vcf info
get_depth <- function(x){
  zero <- 0
  if(substr(as.character(x),1,1)=="."){
    return(zero)
  } else {
    return(as.numeric(unlist(strsplit(as.character(x), ":"))[3]))
  }
}

# import argparse library
suppressMessages(library(argparse))
library(readr)

# Set up cmd line arguments
parser <- ArgumentParser(description="Filter VCF file.")
parser$add_argument("--vcf", action="store", help="Name of VCF file.")
parser$add_argument("--out", action="store", help="Name of output VCF file.")
parser$add_argument("--minQ", action="store", help="min QUAL score for a site.", default="100")
parser$add_argument("--minDP", action="store", help="min read depth for gentype call.", default="5")
parser$add_argument("--missing", action="store", help="max number of missing individuals to include site.")

# Assign parsed arguments to local variables
args    <- parser$parse_args()
vcf_in  <- args$vcf
out     <- args$out
minQ    <- as.numeric(args$minQ)
minDP   <- as.numeric(args$minDP)
missing <- as.numeric(args$missing)

cat("\nReading in VCF file...\n")
vcf <- read.table(vcf_in, header=F, stringsAsFactors = F)

cat("Filtering non-biallelic sites...\n")
biallelic <- apply(as.matrix(vcf[,5]), 1, function(x) nchar(x) == 1)
cat("  Number of non-biallelic sites: ", sum(!biallelic), "\n", sep='')

cat("Filtering sites based on QUAL score...\n")
quality <- (vcf[,6] > minQ)
cat("  Number of poor QUAL sites: ", sum(!quality), "\n", sep='')

cat("Filtering genotypes based on depth...\n")
bad_genos <- 0
for(j in 10:ncol(vcf)){
  idx <- (apply(as.matrix(vcf[,j]), 1, function(x) get_depth(x)) < minDP)
  vcf[idx,j] <- "./."
  bad_genos <- bad_genos + sum(idx)
}
cat("  Number of low depth genotypes: ", bad_genos, "\n", sep='')

cat("Filtering sites based on missing data...\n")
present <- apply(as.matrix(vcf[,10:ncol(vcf)]), 1, function(x) sum(check_missing(x)) < missing)
cat("  Number of sites with greater than ", missing, " missing individuals: ", sum(!present), "\n", sep='')

keepers <- biallelic & quality & present

cat("Writing output VCF...\n")
write_tsv(vcf[keepers,], out, col_names=F)

cat("\n****\n")
cat("Original number of variants:        ", nrow(vcf), "\n", sep='')
cat("Number of variants after filtering: ", sum(keepers), "\n", sep='')
cat("****\n\n")
cat("Done.\n\n")

end <- proc.time() - start
cat("Elapsed time: ", end[3], " seconds.\n\n", sep='')
