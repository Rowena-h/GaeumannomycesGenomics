library(tidyverse)
library(topGO)
library(GO.db)

#Read in orthogroup data
load("R:/GaeumannomycesGenomics/07_comparative_genomics/orthogroup-matrices-2023-07-25.RData")

#Read in strain metadata
metadata <- read.csv(paste0("R:/GaeumannomycesGenomics/05_phylogenomics/raxmlng/metadata.csv"))

#Read in AHRD annotations
ahrds <- do.call("rbind", lapply(
  c(Sys.glob(paste0("S:/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/*/ahrd_output.csv"))),
  function(fn) 
    data.frame(read.csv(fn, sep="\t", header=TRUE, quote="",
                        fill=TRUE, comment.char="#"))
))

#Summarise copy-number and filter out TE genes
copynum.df <- data.frame(strain=colnames(orthogroups.count),
                         as.data.frame(t(orthogroups.count))) %>%
  gather("hog", "num", -strain) %>%
  mutate(biotype=orthogroups.stats$biotype[
    match(hog, orthogroups.stats$orthogroup)
    ]) %>%
  filter(is.na(biotype))

#Filter for high copy-number (HCN) genes
hcn.df <- copynum.df %>%
  filter(num > 10) %>%
  mutate(clade=metadata$clade[match(strain, metadata$strain)])

#Get gene IDs and GO terms for all genes
all.genes <- orthogroups %>%
  #filter(!if_any(everything(.), ~ . == "")) %>%
  rownames_to_column("hog") %>%
  mutate(biotype=orthogroups.stats$biotype[
    match(hog, orthogroups.stats$orthogroup)
    ]) %>%
  filter(is.na(biotype)) %>%
  dplyr::select(-biotype) %>%
  gather(strain, gene, -hog) %>%
  separate_rows(sep=", ", gene) %>%
  mutate(go=ahrds$Gene.Ontology.Term[
    match(gene, sub("-", "", ahrds$Protein.Accession))
    ]) %>%
  filter(go != "")

#Get gene IDs GO terms for HCN genes
hcn.genes <- hcn.df %>%
  left_join(all.genes, by=c("hog", "strain"))

#For each strain
for (strain in unique(hcn.genes$strain)) {
  
  universe <- all.genes[all.genes$strain == strain,]
  test.set <- hcn.genes[hcn.genes$strain == strain,]
  
  #Assign genes of interest amongst gene universe
  geneList <- factor(as.integer(universe$gene %in% test.set$gene))
  names(geneList) <- universe$gene
  
  #Format for topGO
  geneID2GO <- setNames(as.list(as.character(universe$go)), nm=universe$gene)
  
  #Run GO enrichment
  myGOdata <- new("topGOdata", description="Copy-number", ontology="BP",
                  allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
  
  #Fisher's exact test for significance
  resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
  
  #Summarise significant results
  sig.df <- GenTable(myGOdata, raw.p.value=resultFisher, topNodes=length(resultFisher@score)) %>%
    filter(raw.p.value < 0.05) %>%
    mutate(strain=strain)
  
  assign(paste0("sig.df.", strain), sig.df)
  
}

#Combine results
bind_rows(mget(ls(pattern="sig.df.")))
