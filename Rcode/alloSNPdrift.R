Rcpp::sourceCpp("sim_read_data.cpp")
Rcpp::sourceCpp("alloSNPdrift.cpp")

simAlloSNP <- function(f1, ploidy1, ploidy2, nInd, nLoci, coverage, errorDenom=200){

  ploidy = ploidy1 + ploidy2

  freqs2 <- runif(nLoci, 0.05, 0.95)

  genos1 <- matrix(apply(as.matrix(f1), 1, function(x) rbinom(nInd, ploidy1, x)), nrow=nInd, ncol=nLoci)
  genos2 <- matrix(apply(as.matrix(freqs2), 1, function(x) rbinom(nInd, ploidy2, x)), nrow=nInd, ncol=nLoci)
  genos <- genos1 + genos2
  #error <- rbeta(nLoci, 1, errorDenom)
  #totReads <- matrix(rpois(nInd*nLoci, coverage),nrow=nInd, ncol=nLoci)
  #refReads <- simRefReads(totReads, genos, ploidy, error)
  res <- simReads(genos, coverage, ploidy)


  return(list("freqs2"=freqs2,
              "genos1"=genos1,
              "genos2"=genos2,
              "totReads"=res$totalReads,
              "refReads"=res$referenceReads,
              "error"=res$error))

}

# Read in parameters passed at command line
args <- commandArgs(trailingOnly = TRUE)

# Read in ploidy levels for subgenomes
p1 <- as.integer(args[1])
p2 <- as.integer(args[2])

# Number of generations
gen <- as.integer(args[3])

# Number of individuals
N <- as.integer(args[4])

# Sequencing coverage
cover <- as.integer(args[5])

# Outfile prefix
outprefix <- as.character(args[6])

freqs_file <- as.character(args[7])

freqs <- as.matrix(read.table(freqs_file))
f <- freqs


#####################
# Constant pop size #
#####################

fConst <- simConst(f, gen)

# Simulate data
cat("Simulating const pop size...\t")
datConstant <- simAlloSNP(fConst, p1, p2, N, 10000, cover)
cat("Done\n")

# Write everything to file
write.table(fConst, file=paste(outprefix, "-const-freqs1.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datConstant$freqs2, file=paste(outprefix, "-const-freqs2.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datConstant$genos1, file=paste(outprefix, "-const-genos1.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datConstant$genos2, file=paste(outprefix, "-const-genos2.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datConstant$totReads, file=paste(outprefix, "-const-tot.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datConstant$refReads, file=paste(outprefix, "-const-ref.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datConstant$error, file=paste(outprefix, "-const-err.txt", sep=""), quote=F, row.names=F, col.names=F)

############################
# Now sims with pop growth #
############################

fGrowth <- simGrowth(f, gen)

# Simulate data
cat("Simulating with pop growth...\t")
datGrowth <- simAlloSNP(fGrowth, p1, p2, N, 10000, cover)
cat("Done\n")

# Write everything to file
write.table(fGrowth, file=paste(outprefix, "-growth-freqs1.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datGrowth$freqs2, file=paste(outprefix, "-growth-freqs2.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datGrowth$genos1, file=paste(outprefix, "-growth-genos1.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datGrowth$genos2, file=paste(outprefix, "-growth-genos2.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datGrowth$totReads, file=paste(outprefix, "-growth-tot.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datGrowth$refReads, file=paste(outprefix, "-growth-ref.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(datGrowth$error, file=paste(outprefix, "-growth-err.txt", sep=""), quote=F, row.names=F, col.names=F)

warnings()