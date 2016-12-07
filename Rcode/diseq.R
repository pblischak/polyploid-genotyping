# Source cpp files for simulation functions
#require(Rcpp)
Rcpp::sourceCpp("diseq.cpp")
Rcpp::sourceCpp("alloSNP.cpp")
Rcpp::sourceCpp("sim_read_data.cpp")

# Load function for simulating from the diseq model
simDiseq <- function(freqs, diseq, ploidy, nInd, coverage, errorDenom=200){

  #error <- rbeta(length(freqs), 1, errorDenom)
  genotypes <- simDiseqGenos(freqs, diseq, nInd, ploidy);
  #totReads <- matrix(rpois(nInd*length(freqs), coverage),nrow=nInd, ncol=length(freqs))
  
  simData <- simReads(genotypes, coverage, ploidy);

  missing <- simData$totalReads==0
  simData$totalReads[missing] <- -9
  simData$referenceReads[missing] <- -9
  genotypes[missing] <- -9

  return(list("genos"=genotypes, "totReads"=simData$totalReads,
              "refReads"=simData$referenceReads, "error"=simData$error))

}

# Read in parameters passed at command line
args <- commandArgs(trailingOnly = TRUE)

# per-individual inbreedin coeff
F_i <- as.numeric(args[1])

# ploidy level
m_i <- as.integer(args[2])

# number of individuals
N <- as.integer(args[3])

# sequencing coverage
cover <- as.integer(args[4])

# prefix for the outfile
outprefix <- as.character(args[5])

# simulate allele frequencies and write to file
p_ell <- runif(10000, 0.05, 0.95)
write.table(p_ell, file=paste(outprefix, "-freqs.txt", sep=""), quote=F, row.names=F, col.names=F)

# sim data and write to file
dat <- simDiseq(p_ell, rep(F_i, N), m_i, N, cover)
write.table(dat$genos, file=paste(outprefix, "-genos.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(dat$totReads, file=paste(outprefix, "-tot.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(dat$refReads, file=paste(outprefix, "-ref.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(dat$error, file=paste(outprefix, "-err.txt", sep=""), quote=F, row.names=F, col.names=F)