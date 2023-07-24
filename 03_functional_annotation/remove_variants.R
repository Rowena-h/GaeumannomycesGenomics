#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

#Test if there are three arguments: if not, return an error
if (length(args)<3) {
  stop("Three arguments must be supplied: a gff3 file, and output prefix and the directory to save the edited gff3 to.", call.=FALSE)
} 

gff <- read.delim(args[1], header=F, comment.char="#")

#Extract IDs
gff$V10 <- ifelse(
  gff$V3 != "gene",
  sub(".*\\.", "",
      sub("^([^\\.]*\\.[^\\.]*).*", "\\1",
          sub(";.*", "", gff$V9))),
  NA
)

#Filter for the first isoform
gff <- gff[-which(gff$V10 > 1),]

#Write to file
cat("##gff-version 3\n", file=paste0(args[3], args[2], ".gff3"))
write.table(
  gff[,1:9],
  paste0(args[3], "/", args[2], ".gff3"),
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE
)
