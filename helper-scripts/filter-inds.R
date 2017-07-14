#!/usr/bin/Rscript

# Return a logical vector specifying if a column is to be included (TRUE)
# or excluded (FALSE)
get_missing <- function(col, cutoff){
  if(mean(is.na(col)) > cutoff){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

filter <- function(tot_file, alt_file, missing_cutoff = 0.5,
                   transpose = TRUE, missing_string = "-9"){
  tot <- as.matrix(read.table(tot_file, na.strings = missing_string))
  alt <- as.matrix(read.table(alt_file, na.strings = missing_string))
  filter <- apply(tot, 2, get_missing, missing_cutoff)

  cat("Number of individuals: ", dim(tot[,filter])[2], "\n")
  cat("Number of loci: ", dim(tot[,filter])[1], "\n")

  if(transpose){
    write.table(t(tot[,filter]), file=paste("filtered-", tot_file, sep=""),
              quote=FALSE, sep="\t", na="-9", row.names=FALSE, col.names=FALSE)
    write.table(t(alt[,filter]), file=paste("filtered-", alt_file, sep=""),
              quote=FALSE, sep="\t", na="-9", row.names=FALSE, col.names=FALSE)
  } else {
    write.table(tot[,filter], file=paste("filtered-", tot_file, sep=""),
              quote=FALSE, sep="\t", na="-9", row.names=FALSE, col.names=FALSE)
    write.table(alt[,filter], file=paste("filtered-", alt_file, sep=""),
              quote=FALSE, sep="\t", na="-9", row.names=FALSE, col.names=FALSE)
  }
}

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 2){
  stop("No input file names specified.")
} else {
  filter(args[1], args[2], as.numeric(args[3]), as.logical(args[4]), args[5])
}
