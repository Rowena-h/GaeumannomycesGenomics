#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)

#Test if there are three arguments: if not, return an error
if (length(args)<3) {
  stop("Three arguments must be supplied: a gff3 file, and output prefix and the directory to save the edited gff3 to.", call.=FALSE)
}

library(tidyverse)

gff <- read.delim(args[1], header=F, comment.char="#")

#Extract IDs and remove existing notes
gff.id <- gff %>%
  mutate(ID=sub("\\..*", "", sub(";.*", "", V9)),
         V9=sub("Note=.*;", "", V9))

#Get genes with variants
variants <- gff.id %>%
  filter(V3 == "mRNA") %>%
  group_by(ID) %>%
  summarise(num.variants=n()) %>%
  filter(num.variants > 1)

#Add note for variants
gff.variant.note <- gff.id %>%
  mutate(V9=ifelse(
    ID %in% variants$ID & V3 != "gene",
    sub("$", ";Note=alternatively spliced", V9),
    V9
  ))

#Write to file
cat("##gff-version 3\n", file=paste0(args[3], args[2], ".gff3"))
write.table(
  gff.variant.note[,1:9],
  paste0(args[3], "/", args[2], ".gff3"),
  sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE
)

