setwd("R://GaeumannomycesGenomics/04_phylogenetic_classification/")

library(ape)
library(tidyverse)
library(ggtree)
library(ggpubr)

##gdo tree

#Read in tree
gdo.tree <- read.tree("raxmlng/gaeumannomyces/gaeumannomyces_gdo.raxml.support")
#Root
gdo.tree <- root(gdo.tree,
                      c("Gaeumannomyces_graminis_var._tritici_haplotype_17", "Gt-3aA1", "Gt-CB1"),
                      edgelabel=TRUE, resolve.root=TRUE)
#Truncate branch
gdo.shortened.edge <- gdo.tree$edge[which.max(gdo.tree$edge.length), 2]
gdo.tree$edge.length[which.max(gdo.tree$edge.length)] <- gdo.tree$edge.length[which.max(gdo.tree$edge.length)] / 3

#Get metadata
gaeu.metadata <- read.csv("raxmlng/gaeumannomyces/metadata.csv")
#Format tip labels
gaeu.metadata$new.label <- ifelse(gaeu.metadata$own == "Y",
                                 paste0('paste(bold("', gaeu.metadata$name, '"))'),
                                 paste0('paste("', gaeu.metadata$name, '")'))

#Plot tree
gg.gdo.tree <- ggtree(gdo.tree, linetype=NA) %<+% gaeu.metadata +
  geom_tree(aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
            show.legend=FALSE) +
  xlim(0, 0.045) +
  scale_colour_manual(values="grey",
                      na.value="black") +
  geom_label2(aes(x=branch, subset=node == gdo.shortened.edge),
              label="//",
              size=3,
              label.padding=unit(0, "pt"),
              label.size=0) +
  geom_tiplab(aes(label=new.label),
              size=2,
              parse=TRUE)

gg.gdo.tree <- gg.gdo.tree +
  annotate(geom="label",
           label="gdo",
           x=max(na.omit(gg.gdo.tree$data$x))*0.3,
           y=max(na.omit(gg.gdo.tree$data$y))*0.92,
           size=3)

#Make dataframe of cluster nodes
gdo.clades.df <- data.frame(
  clade=unique(gaeu.metadata %>%
                 filter(clade != "") %>%
                 pull(clade)),
  node=NA
)

#Find the most recent common ancestor for each clade
for (j in 1:length(gdo.clades.df$clade)) {
  
  gdo.clades.df$node[j] <- MRCA(gdo.tree,
                                gaeu.metadata$tip[gaeu.metadata$clade == gdo.clades.df$clade[j]])
  
}

#Add clade labels
gg.gdo.tree <- gg.gdo.tree +
  geom_cladelab(data=gdo.clades.df,
                mapping=aes(node=node, label=clade),
                fontsize=3,
                barsize=0.8,
                offset=0.012,
                offset.text=0.001)



##ITS2 tree

#Read in tree
its2.tree <- read.tree("raxmlng/gaeumannomyces/gaeumannomyces_ITS2.raxml.support")
#Root
its2.tree <- root(its2.tree,
                 c("Gaeumannomyces_avenae", "Gt-3aA1", "Gt-CB1"),
                 edgelabel=TRUE, resolve.root=TRUE)

#Plot tree
gg.its2.tree <- ggtree(its2.tree, linetype=NA) %<+% gaeu.metadata +
  geom_tree(aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
            show.legend=FALSE) +
  scale_x_reverse(limits=c(0.07, 0)) +
  scale_colour_manual(values="grey",
                      na.value="black") +
  geom_tiplab(aes(label=new.label),
              size=2,
              hjust=1,
              offset=-0.001,
              parse=TRUE)

#Add clade labels
gg.its2.tree <- gg.its2.tree +
  geom_strip("Gt-4e", "Gaeumannomyces_tritici_B",
             label="B",
             fontsize=3,
             barsize=0.8,
             offset=0.01,
             offset.text=0.005) +
  geom_strip("Gt-19d1", "Gt-8d",
             label="A",
             fontsize=3,
             barsize=0.8,
             offset=-0.01,
             offset.text=0.005) +
  geom_strip("Gt-CB1", "Gt-3aA1",
             label="avenae",
             fontsize=3,
             barsize=0.8,
             offset=-0.03,
             offset.text=0.02) +
  annotate(geom="label",
           label="ITS2",
           x=max(na.omit(gg.its2.tree$data$x))*0.3,
           y=max(na.omit(gg.its2.tree$data$y))*0.92,
           size=3)

#Write to file
tiff(file=paste0("gaeumannomyces_groups-", Sys.Date(), ".tiff"),
     height=4, width=6, units="in",
     res=600, compression="lzw")
ggarrange(gg.gdo.tree, gg.its2.tree)
dev.off()

