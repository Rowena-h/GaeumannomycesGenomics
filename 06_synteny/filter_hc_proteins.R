#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

#Test if there are three arguments: if not, return an error
if (length(args)<3) {
  stop("Three arguments must be supplied: a gff3 file, an output prefix and the directory to save the edited gff3 to.", call.=FALSE)
}

library(tidyverse)

gff <- read.delim(args[1], header=F, comment.char="#")

#Extract IDs, biotypes and confidence
gff.id <- gff %>%
  mutate(ID=sub("\\..*", "", sub(";.*", "", V9)),
         biotype=ifelse(V3 == "gene", sub(";.*", "", sub(".*biotype=", "", V9)), NA),
         confidence=ifelse(V3 == "gene", sub(";.*", "", sub(".*confidence=", "", V9)), NA))

#Get high confidence protein coding genes
hc.protein.ids <- gff.id %>%
  filter(biotype == "protein_coding_gene",
         confidence == "High") %>%
  pull(ID)

#Get genes with variants
gff.hc.proteins <- gff.id %>%
  filter(ID %in% hc.protein.ids)

cat(paste0(strain, ": ", length(unique(gff.hc.proteins$ID)), " high confidence protein coding gene models (out of ", length(unique(gff.id$ID)), ")"))

#Write to file
cat("##gff-version 3\n", file=paste0(args[3], "/", args[2], ".gff3"))
write.table(
  gff.hc.proteins[,1:9],
  paste0(args[3], "/", args[2], ".gff3"),
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE
)
