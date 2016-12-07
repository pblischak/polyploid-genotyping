Rcpp::sourceCpp("sim_read_data.cpp")
#Rcpp::sourceCpp("alloSNP.cpp")

simAlloSNP <- function(ploidy1, ploidy2, nInd, nLoci, coverage, errorDenom=200){

  ploidy = ploidy1 + ploidy2

  freqs1 <- runif(nLoci, 0.05, 0.95)
  freqs2 <- runif(nLoci, 0.05, 0.95)

  genos1 <- matrix(apply(as.matrix(freqs1), 1, function(x) rbinom(nInd, ploidy1, x)), nrow=nInd, ncol=nLoci)
  genos2 <- matrix(apply(as.matrix(freqs2), 1, function(x) rbinom(nInd, ploidy2, x)), nrow=nInd, ncol=nLoci)
  genos <- genos1 + genos2
  #error <- rbeta(nLoci, 1, errorDenom)
  #totReads <- matrix(rpois(nInd*nLoci, coverage),nrow=nInd, ncol=nLoci)
  #refReads <- simRefReads(totReads, genos, ploidy, error)
  res <- simReads(genos, coverage, ploidy)


  return(list("freqs1"=freqs1,
              "freqs2"=freqs2,
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

# Number of individuals
N <- as.integer(args[3])

# Sequencing coverage
cover <- as.integer(args[4])

# Outfile prefix
outprefix <- as.character(args[5])

# Simulate data
dat <- simAlloSNP(p1, p2, N, 10000, cover)

# Write everything to file
write.table(dat$freqs1, file=paste(outprefix, "-freqs1.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(dat$freqs2, file=paste(outprefix, "-freqs2.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(dat$genos1, file=paste(outprefix, "-genos1.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(dat$genos2, file=paste(outprefix, "-genos2.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(dat$totReads, file=paste(outprefix, "-tot.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(dat$refReads, file=paste(outprefix, "-ref.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(dat$error, file=paste(outprefix, "-err.txt", sep=""), quote=F, row.names=F, col.names=F)
