##Script to produce input files for lifestyle test:
##https://github.com/fantin-mesny/Effect-Of-Biological-Categories-On-Genomes-Composition

library(ape)
library(tidyverse)

#Directory paths
dir.phylo <- "R:/GaeumannomycesGenomics/05_phylogenomics/"
dir.comp <- "R:/GaeumannomycesGenomics/07_comparative_genomics/"


#Read in metadata
metadata <- read.csv(paste0(dir.phylo, "raxmlng/metadata.csv"))


## Format species tree ##

#Read in tree
tree <- read.tree(paste0(dir.phylo, "raxmlng/gaeumannomyces_concat.raxml.support"))
tree$tip.label <- sub("_EIv1", "", tree$tip.label)

#Remove outgroup from tree and write to file
outgroups <- c("GCA_000193285.1_Mag_poae_ATCC_64411_V1_protein.faa")
tree <- root(tree, outgroups, edgelabel=TRUE, resolve.root=TRUE)
tree.ingroup <- drop.tip(
  tree, c(outgroups, "GCF_000145635.1_Gae_graminis_V2_protein.faa_XP")
)
write.tree(tree.ingroup,
           paste0(dir.comp, "lifestyle_permanova/species_tree_ingroup.tre"))


## Format gene abundance matrices ##

#For each gene type...
for (gene in c("orthogroup", "CAZyme", "CSEP")) {
  
  #Read in abundance matrix
  matrix <- read.csv(paste0(dir.comp, "lifestyle_permanova/", gene, "-count.csv")) %>%
    column_to_rownames(var="X") %>%
    filter(rowSums(across())>0)
  
  colnames(matrix) <- sub("\\.", "-", colnames(matrix))
  
  #Transpose
  matrix <- as.data.frame(t(matrix))
  
  #Add column with names
  matrix$genome <- rownames(matrix)
  
  #Add column with lifestyle
  matrix$lifestyle <- ifelse(grepl("Gh", matrix$genome),
                             "nonpathogen", "pathogen")
  #Reorder columns
  matrix <- matrix %>%
    select(genome, lifestyle, everything())
  
  #Write to file
  write.csv(matrix,
            paste0(dir.comp, "lifestyle_permanova/", gene, "_abundance_matrix.csv"),
            row.names=FALSE, quote=FALSE)
  
}


## Format BGC abundance matrix ##

bgc.dir <- "S:/functional_annotation/002.5_big-scape/gaeumannomyces/"

#Read in BGC data
bgc.df <- read.csv(
  Sys.glob(paste0(bgc.dir, "network_files/*/Network_Annotations_Full.tsv")),
  sep="\t", header=TRUE
)

bgc.clusters.df <- do.call("rbind", lapply(
  c(Sys.glob(paste0(bgc.dir, "network_files/*/PKSI/PKSI_clustering_c0.30.tsv")),
    Sys.glob(paste0(bgc.dir, "network_files/*/NRPS/NRPS_clustering_c0.30.tsv")),
    Sys.glob(paste0(bgc.dir, "network_files/*/PKS-NRP_Hybrids/PKS-NRP_Hybrids_clustering_c0.30.tsv")),
    Sys.glob(paste0(bgc.dir, "network_files/*/PKSother/PKSother_clustering_c0.30.tsv")),
    Sys.glob(paste0(bgc.dir, "network_files/*/Terpene/Terpene_clustering_c0.30.tsv"))),
  function(fn) 
    data.frame(read.csv(fn, sep="\t", header=TRUE))
))

#Make abundance matrix
bgc.matrix <- bgc.clusters.df %>%
  rename(BGC="X.BGC.Name") %>%
  inner_join(bgc.df) %>%
  mutate(taxon=sub("_EI_v1.1.*", "", Accession.ID)) %>%
  group_by(taxon, Family.Number) %>%
  summarise(num=n()) %>%
  pivot_wider(names_from=Family.Number,
              values_from=num,
              values_fill=0) %>%
  rename(genome="taxon")

#Add column with lifestyle
bgc.matrix$lifestyle <- ifelse(grepl("Gh", bgc.matrix$genome),
                               "nonpathogen", "pathogen")
#Reorder columns
bgc.matrix <- bgc.matrix %>%
  select(genome, lifestyle, everything())

#Write to file
write.csv(bgc.matrix, 
          paste0(dir.comp, "lifestyle_permanova/BGC_abundance_matrix.csv"),
          row.names=FALSE, quote=FALSE)

#Submit test
system("sbatch permanova.sh")