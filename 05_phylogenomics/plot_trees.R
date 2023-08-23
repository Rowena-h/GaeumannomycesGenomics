library(ape)
library(tidyverse)
library(ggtree)
library(aplot)
library(ggplotify)
library(patchwork)
library(scales)

## gdo tree ##

#Read in tree
gdo.tree <- read.tree("R://GaeumannomycesGenomics/04_phylogenetic_classification/raxmlng/gaeumannomyces_gdo.raxml.support")
#Root
gdo.tree <- root(gdo.tree,
                 c("Gaeumannomyces_graminis_var._tritici_haplotype_17", "Gt-3aA1", "Gt-CB1"),
                 edgelabel=TRUE, resolve.root=TRUE)
#Truncate branch
gdo.shortened.edge <- gdo.tree$edge[which.max(gdo.tree$edge.length), 2]
gdo.tree$edge.length[which.max(gdo.tree$edge.length)] <- gdo.tree$edge.length[which.max(gdo.tree$edge.length)] / 3

#Get metadata
markers.metadata <- read.csv("R://GaeumannomycesGenomics/04_phylogenetic_classification/raxmlng/metadata.csv")
#Format tip labels
markers.metadata$new.label <- ifelse(
  markers.metadata$own == "Y",
  paste0('paste(bold("', markers.metadata$name, '"))'),
  paste0('paste("', markers.metadata$name, '")')
)

#Plot gdo tree
gg.gdo.tree <- ggtree(gdo.tree, linetype=NA) %<+% markers.metadata +
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
              offset=0.001,
              parse=TRUE)

gg.gdo.tree <- gg.gdo.tree +
  annotate(geom="label",
           label="gdo",
           x=max(na.omit(gg.gdo.tree$data$x))*0.3,
           y=max(na.omit(gg.gdo.tree$data$y))*0.92,
           size=3)

#Make dataframe of clade nodes
gdo.clades.df <- data.frame(
  clade=unique(markers.metadata %>%
                 filter(clade != "") %>%
                 pull(clade)),
  node=NA
)

#Find the most recent common ancestor for each clade
for (j in 1:length(gdo.clades.df$clade)) {
  
  gdo.clades.df$node[j] <- MRCA(gdo.tree,
                                markers.metadata$tip[markers.metadata$clade == gdo.clades.df$clade[j]])
  
}

#Add clade labels
gg.gdo.tree <- gg.gdo.tree +
  geom_cladelab(data=gdo.clades.df,
                mapping=aes(node=node, label=clade),
                fontsize=3,
                barsize=0.8,
                offset=0.012,
                offset.text=0.001)


## ITS2 tree ##

#Read in tree
its2.tree <- read.tree("R://GaeumannomycesGenomics/04_phylogenetic_classification/raxmlng/gaeumannomyces_ITS2.raxml.support")
#Root
its2.tree <- root(its2.tree,
                  c("Gaeumannomyces_avenae", "Gt-3aA1", "Gt-CB1"),
                  edgelabel=TRUE, resolve.root=TRUE)

#Plot ITS2 tree
gg.its2.tree <- ggtree(its2.tree, linetype=NA) %<+% markers.metadata +
  geom_tree(aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
            show.legend=FALSE) +
  scale_x_reverse(limits=c(0.07, 0)) +
  scale_colour_manual(values="grey",
                      na.value="black") +
  geom_tiplab(aes(label=new.label),
              size=2,
              hjust=1,
              offset=-0.001,
              parse=TRUE) +
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
             offset.text=0.02)

gg.its2.tree <- gg.its2.tree +
  annotate(geom="label",
           label="ITS2",
           x=max(na.omit(gg.its2.tree$data$x))*0.3,
           y=max(na.omit(gg.its2.tree$data$y))*0.92,
           size=3)


## SPECIES TREE

#Read in tree
species.tree <- read.tree(paste0("R://GaeumannomycesGenomics/05_phylogenomics/raxmlng/gaeumannomyces_concat.raxml.support"))

#Root tree
species.tree <- root(species.tree, "GCA_000193285.1_Mag_poae_ATCC_64411_V1_protein.faa",
                     edgelabel=TRUE, resolve.root=TRUE)

#Truncate excessively long branch
shortened.edge <- species.tree$edge[which.max(species.tree$edge.length), 2]
species.tree$edge.length[which.max(species.tree$edge.length)] <- species.tree$edge.length[which.max(species.tree$edge.length)] / 3

#Read in tree metadata
metadata <- read.csv(paste0("R://GaeumannomycesGenomics/05_phylogenomics/raxmlng/metadata.csv"))
#Format tip labels
metadata$new.label <- ifelse(
  metadata$own == "Y",
  paste0('paste(bolditalic("', metadata$name, '"), bold(" ', metadata$new.strain, '"))'),
  paste0('paste(italic("', metadata$name, '"), " ', metadata$strain, '")')
)

#Plot species tree
gg.species.tree <- ggtree(species.tree, linetype=NA) %<+% metadata

#Make dataframe of Gt type nodes
types.df <- data.frame(
  type=unique(metadata %>%
                filter(type != "") %>%
                pull(type)),
  node=NA
)

#Find the most recent common ancestor for each type
for (i in 1:length(types.df$type)) {
  
  types.df$node[i] <- 
    MRCA(species.tree,
         metadata$tip[metadata$type == types.df$type[i]])
  
}

#Make dataframe of clade nodes
clades.df <- data.frame(clade=unique(metadata$clade[which(metadata$clade != "outgroup")]),
                        node=NA)

#Find the most recent common ancestor for each clade
for (i in 1:length(clades.df$clade)) {
  
  clades.df$node[i] <- 
    MRCA(gg.species.tree,
         metadata$tip[metadata$clade == clades.df$clade[i]])
  
}

#Make alternated coding for highlights on tree
clades.df <- clades.df[match(
  clades.df$clade,
  unique(gg.species.tree$data %>%
           filter(clade != "outgroup" & !is.na(clade) & isTip == "TRUE") %>%
           arrange(y) %>%
           pull(clade))
),]

clades.df$box <- rep(c(0,1), length.out=length(clades.df$clade))

#Add clade labels and highlights
gg.species.tree <- gg.species.tree +
  geom_highlight(data=clades.df, 
                 aes(node=node, fill=as.factor(box)),
                 alpha=1,
                 align="right",
                 extend=0.035,
                 show.legend=FALSE) +
  geom_cladelab(data=types.df,
                mapping=aes(node=node, label=type),
                fontsize=3,
                barsize=0.8,
                align=TRUE,
                offset=0.025,
                offset.text=0.001) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  geom_tree(aes(linetype=ifelse(as.numeric(label) < 70, "insig", NA)),
            show.legend=FALSE) +
  xlim(0, 0.1) +
  scale_linetype_manual(values="11", 
                        na.value="solid") +
  geom_tiplab(aes(label=new.label),
              offset=0.001,
              size=3,
              parse=TRUE) +
  annotate(geom="label",
           label="2,359 loci",
           x=max(na.omit(gg.species.tree$data$x))*0.3,
           y=max(na.omit(gg.species.tree$data$y))*0.92,
           size=3) +
  geom_label2(aes(x=branch, subset=node == shortened.edge),
              label="//",
              label.padding=unit(0, "pt"),
              label.size=0) +
coord_cartesian(clip="off") +
  theme(legend.title=element_blank(),
        legend.position=c(0.25, 0.75),
        legend.margin=margin(0, 0, 0, 0))

#Write to file
pdf(file=paste0("R://GaeumannomycesGenomics/05_phylogenomics/trees-", Sys.Date(), ".pdf"),
    height=7, width=7)
wrap_elements(gg.gdo.tree + gg.its2.tree + 
                  plot_annotation(title="a") &
                  theme(plot.title=element_text(face="bold"))) /
  wrap_elements(gg.species.tree + 
                  plot_annotation(title="b") &
                  theme(plot.title=element_text(face="bold")))
dev.off()


## R3-111-a comparison ##

#Plot bargraph of number of genes
gg.gene.numbers <- ggplot(metadata, aes(y=tip, x=genes)) +
  geom_bar(stat="identity", width=0.6) +
  geom_label(aes(label=comma(genes)), 
             position=position_stack(vjust=0.5),
             size=1.5,
             label.padding=unit(0.1, "lines"),
             colour="#595959") +
  labs(x="Number of gene models") +
  scale_x_continuous(expand=c(0, 0),
                     position="top",
                     labels=comma) +
  theme_minimal() +
  theme(axis.title.x=element_text(size=7, face="bold"),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank())

#Read in number of functionally annotated genes
df.ahrd.annotated <- read.csv("R:/GaeumannomycesGenomics/03_functional_annotation/ahrd_comparison.tsv", sep="\t") %>%
  mutate(new.strain=c("Gh-1B17", "Gh-2C17", "Gt-LH10", "Gt-19d1", "Gt-23d",
                      "Ga-3aA1", "Gt-4e", "Gt-8d", "Ga-CB1", "R3-111a-1"),
         tip=metadata$tip[match(new.strain, metadata$new.strain)],
         proportion.primary=annotated_primary_transcripts/primary_transcripts)

#Plot bargraph of number of functionally annotated genes
gg.annotation.proportion <- ggplot(df.ahrd.annotated, aes(x=proportion.primary, y=tip)) +
  geom_bar(stat="identity", width=0.6) +
  geom_label(aes(label=round(proportion.primary, 2)), 
             position=position_stack(vjust=0.5),
             size=1.5,
             label.padding=unit(0.1, "lines"),
             colour="#595959") +
  labs(x="Proportion of functionally\nannotated gene models") +
  scale_x_continuous(expand=c(0, 0),
                     position="top") +
  theme_minimal() +
  theme(axis.title.x=element_text(size=7, face="bold"),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank())

#Read in OrthoFinder unassigned genes
df.orthofinder.unassigned <- read.csv("S:/orthofinder/gaeumannomyces/Results_Apr28/Orthogroups/Orthogroups_UnassignedGenes.tsv", sep="\t") %>%
  select(-c("Orthogroup", "GCA_000193285.1_Mag_poae_ATCC_64411_V1_protein")) %>%
  gather(strain, unassigned.genes) %>%
  filter(unassigned.genes != "") %>%
  group_by(strain) %>%
  summarise(num.unassigned.genes=n()) %>%
  mutate(new.strain=c("R3-111a-1", "Gh-1B17", "Gh-2C17", "Gt-19d1", "Gt-23d",
                      "Ga-3aA1", "Gt-4e", "Gt-8d", "Ga-CB1", "Gt-LH10"),
         tip=metadata$tip[match(new.strain, metadata$new.strain)])

#Plot bargraph of unassigned genes
gg.unassigned.genes <- ggplot(df.orthofinder.unassigned, aes(y=tip, x=num.unassigned.genes)) +
  geom_bar(stat="identity", width=0.6) +
  geom_label(aes(label=comma(num.unassigned.genes)), 
             position=position_stack(vjust=0.5),
             size=1.5,
             label.padding=unit(0.1, "lines"),
             colour="#595959") +
  scale_x_continuous(expand=c(0, 0),
                     position="top",
                     labels=comma) +
  coord_cartesian(clip="off") +
  labs(x="Number of OrthoFinder\nunassigned gene models") +
  theme_minimal() +
  theme(axis.title.x=element_text(size=7, face="bold"),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank())

#Plot species tree with no branch lengths
gg.species.tree.comp <- ggtree(species.tree, linetype=NA, branch.length="none") %<+% metadata +
  geom_highlight(data=clades.df, 
                 aes(node=node, fill=as.factor(box)),
                 alpha=1,
                 align="right",
                 extend=100,
                 show.legend=FALSE) +
  geom_cladelab(data=types.df,
                mapping=aes(node=node, label=type),
                fontsize=3,
                barsize=0.8,
                align=TRUE,
                offset=7,
                offset.text=0.5) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  geom_tree(aes(linetype=ifelse(as.numeric(label) < 70, "insig", NA)),
            show.legend=FALSE) +
  xlim(0, 15) +
  scale_linetype_manual(values="11", 
                        na.value="solid") +
  geom_tiplab(aes(label=new.label),
              offset=0.1,
              size=3,
              parse=TRUE) +
  geom_label2(aes(x=branch, subset=node == shortened.edge),
              label="//",
              label.padding=unit(0, "pt"),
              label.size=0) +
  coord_cartesian(clip="off") +
  geom_highlight(node=MRCA(species.tree, metadata$tip[metadata$name == "G. tritici" & metadata$own != "Y"]),
                 fill=NA,
                 colour="black",
                 linetype="dashed",
                 alpha=0, extend=28.9,
                 show.legend=FALSE)

#Write to file
pdf(paste0("R://GaeumannomycesGenomics/05_phylogenomics/annotation_comp-", Sys.Date(), ".pdf"), width=8, height=4)

gg.gene.numbers %>% 
  insert_left(gg.species.tree.comp, width=2.7) %>%
  insert_right(gg.annotation.proportion) %>%
  insert_right(gg.unassigned.genes)

dev.off()

