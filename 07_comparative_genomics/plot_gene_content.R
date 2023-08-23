library(ape)
library(tidyverse)
library(aplot)
library(ggtree)
library(scales)
library(vegan)
library(ggpubr)
library(rstatix)
library(matrixStats)
library(eulerr)
library(ggplotify)

## Species tree ##

#Assign outgroup
outgroups <- c("GCA_000193285.1_Mag_poae_ATCC_64411_V1_protein.faa")

#Read in tree metadata
metadata <- read.csv(paste0("R:/GaeumannomycesGenomics/05_phylogenomics/raxmlng/metadata.csv"))
#Format tip labels
metadata$new.label <- ifelse(
  metadata$own == "Y",
  paste0('paste(bolditalic("', metadata$name, '"), bold(" ', metadata$new.strain, '"))'),
  paste0('paste(italic("', metadata$name, '"), " ', metadata$new.strain, '")')
)

#Read in tree
tree <- read.tree("R:/GaeumannomycesGenomics/05_phylogenomics/raxmlng/gaeumannomyces_concat.raxml.support")

#Root tree
tree <- root(tree, outgroups, edgelabel=TRUE, resolve.root=TRUE)

#Remove other Gt tip
tree <- drop.tip(tree, "GCF_000145635.1_Gae_graminis_V2_protein.faa_XP")

#Truncate excessively long branch
shortened.edge <- tree$edge[which.max(tree$edge.length), 2]
tree$edge.length[which.max(tree$edge.length)] <- tree$edge.length[which.max(tree$edge.length)] / 3

#Plot tree
gg.tree <- ggtree(tree, linetype=NA) %<+% metadata

#Make dataframe of Gt type nodes
types.df <- data.frame(
  type=unique(metadata %>%
                filter(type != "") %>%
                pull(type)),
  node=NA
)

#Find the most recent common ancestor for each type
for (j in 1:length(types.df$type)) {
  
  types.df$node[j] <- MRCA(tree,
                           metadata$tip[metadata$type == types.df$type[j] &
                                          metadata$own == "Y"])
  
}

#Make dataframe of clade nodes
clades.df <- data.frame(clade=unique(metadata$clade[which(metadata$clade != "outgroup")]),
                        node=NA)

#Find the most recent common ancestor for each ckade
for (i in 1:length(clades.df$clade)) {
  clades.df$node[i] <- 
    MRCA(gg.tree,
         metadata$tip[metadata$clade == clades.df$clade[i] &
                        metadata$own == "Y"])
}

#Make alternated coding for highlights on tree
clades.df <- clades.df[match(clades.df$clade,
                             unique(gg.tree$data %>%
                                      filter(clade != "outgroup" & !is.na(clade) & isTip == "TRUE") %>%
                                      arrange(y) %>%
                                      pull(clade))),]

clades.df$box <- rep(c(0,1), length.out=length(clades.df$clade))

#Add clade labels and highlights and tip labels
gg.tree <- gg.tree +
  geom_highlight(data=clades.df, 
                 aes(node=node, fill=as.factor(box)),
                 alpha=1, extend=0.2,
                 show.legend=FALSE) +
  geom_cladelab(data=types.df,
                mapping=aes(node=node, label=type),
                fontsize=3,
                barsize=0.8,
                offset=0.032,
                offset.text=0.001) +
  scale_fill_manual(values=c("#F5F5F5", "#ECECEC")) +
  geom_tree(aes(linetype=ifelse(as.numeric(label) < 70, "insig", NA)),
            show.legend=FALSE) +
  scale_linetype_manual(values="11", 
                        na.value="solid") +
  geom_tiplab(aes(label=new.label),
              offset=0.001,
              size=3,
              parse=TRUE) +
  geom_label2(aes(x=branch, subset=node == shortened.edge),
              label="//",
              label.padding=unit(0, "pt"),
              label.size=0) +
  xlim(0, 0.105) +
  coord_cartesian(clip="off") +
  theme(legend.title=element_blank(),
        legend.position=c(0.25, 0.75),
        legend.margin=margin(0, 0, 0, 0))


## Lifestyle gene variance comparison PERMANOVA ##

#For all types of genes...
for (i in c("orthogroup", "CSEP", "CAZyme", "BGC")){
  
  print(i)
  lifestyle.data <- read.csv(
    paste0("R:/GaeumannomycesGenomics/07_comparative_genomics/", i, "/data.csv"),
    row.names="genome"
  )
  #Make distance matrix of orthogroup content
  dist <- vegdist(lifestyle.data, method="jaccard")
  #Read in lifestyle test results
  phy.pca.result <- read.csv(
    paste0("R:/GaeumannomycesGenomics/07_comparative_genomics/", i, "/metadata.csv")
  )
  #Do permanova
  permanova <- adonis2(formula=dist ~ PC1 + PC2 + lifestyle,
                       data=phy.pca.result,
                       permutations=9999)
  
  print(paste0("Phylogeny: ", round(sum(permanova$R2[1:2]) * 100), "%"))
  print(paste0("Lifestyle: ", round(sum(permanova$R2[3]) * 100), "%"))
  
  assign(paste0("permanova.", i), permanova)
  assign(paste0("phy.pca.result.", i), phy.pca.result)
  
}


## Gene content summary ##

# BGCS #

#Read in BiG-SCAPE results
bgc.dir <- "S:/functional_annotation/002.5_big-scape/gaeumannomyces/"
bgc.df <- read.csv(Sys.glob(paste0(bgc.dir, "network_files/*/Network_Annotations_Full.tsv")), sep="\t", header=TRUE)

bgc.clusters.df <- do.call("rbind", lapply(
  c(Sys.glob(paste0(bgc.dir, "network_files/*/PKSI/PKSI_clustering_c0.30.tsv")),
    Sys.glob(paste0(bgc.dir, "network_files/*/NRPS/NRPS_clustering_c0.30.tsv")),
    Sys.glob(paste0(bgc.dir, "network_files/*/PKS-NRP_Hybrids/PKS-NRP_Hybrids_clustering_c0.30.tsv")),
    Sys.glob(paste0(bgc.dir, "network_files/*/PKSother/PKSother_clustering_c0.30.tsv")),
    Sys.glob(paste0(bgc.dir, "network_files/*/Terpene/Terpene_clustering_c0.30.tsv"))),
  function(fn) 
    data.frame(read.csv(fn, sep="\t", header=TRUE))
))

#Summarise abundance of BGC families
bgc.num.df <- bgc.clusters.df %>%
  dplyr::rename(BGC="X.BGC.Name") %>%
  inner_join(bgc.df) %>%
  mutate(strain=sub("_EI_v1.1.*", "", Accession.ID)) %>%
  group_by(strain, Family.Number) %>%
  summarise(num=n())

#Create abundance matrix
bgc.abundance.mat <- bgc.num.df %>%
  pivot_wider(names_from=Family.Number,
              values_from=num,
              values_fill=0)

#Filter for only Gt strains
bgc.abundance.mat.gt <- bgc.abundance.mat %>%
  filter(strain %in% metadata$strain[metadata$name == "G. tritici"]) %>%
  column_to_rownames("strain")

bgc.abundance.mat.gt <- bgc.abundance.mat.gt[,colSums(bgc.abundance.mat.gt) > 0]

bgc.num.df.gt <- bgc.num.df %>%
  filter(strain %in% metadata$strain[metadata$name == "G. tritici"]) %>%
  mutate(category=NA)

#Categorise Gt BGCs as core/accessory/specific
for (i in 1:length(colnames(bgc.abundance.mat.gt))) {
  
  if (length(which(bgc.abundance.mat.gt[,i] == 0)) == (length(rownames(bgc.abundance.mat.gt)) - 1)) {
    bgc.num.df.gt$category[bgc.num.df.gt$Family.Number == colnames(bgc.abundance.mat.gt)[i]] <-
      "specific"
  }
  if (length(which(bgc.abundance.mat.gt[,i] == 0)) == 0) {
    bgc.num.df.gt$category[bgc.num.df.gt$Family.Number == colnames(bgc.abundance.mat.gt)[i]] <-
      "core"
  }
  if (length(which(bgc.abundance.mat.gt[,i] > 0)) == (length(rownames(bgc.abundance.mat.gt)) - 1)) {
    bgc.num.df.gt$category[bgc.num.df.gt$Family.Number == colnames(bgc.abundance.mat.gt)[i]] <-
      "softcore"
  }
  if (length(which(bgc.abundance.mat.gt[,i] == 0)) < (length(rownames(bgc.abundance.mat.gt)) - 1) &&
      length(which(bgc.abundance.mat.gt[,i] == 0)) > 1) {
    bgc.num.df.gt$category[bgc.num.df.gt$Family.Number == colnames(bgc.abundance.mat.gt)[i]] <-
      "accessory"
  }
  
}


# CSEPs/CAZymes #

#Read in orthogroup data
load("R:/GaeumannomycesGenomics/07_comparative_genomics/orthogroup-matrices-2023-07-25.RData")

#Filter for only Gt strains
orthogroups.count.gt <- orthogroups.count %>%
  mutate(orthogroup=rownames(orthogroups.count),
         biotype=orthogroups.stats$biotype[match(orthogroup, orthogroups.stats$orthogroup)]) %>%
  filter(is.na(biotype)) %>%
  select(colnames(orthogroups.count)[colnames(orthogroups.count) %in% metadata$strain[metadata$name == "G. tritici"]]) %>%
  filter(rowSums(across())>0)

CSEP.count.gt <- CSEP.count %>%
  select(colnames(orthogroups.count)[colnames(orthogroups.count) %in% metadata$strain[metadata$name == "G. tritici"]]) %>%
  filter(rowSums(across())>0)

CAZyme.count.gt <- CAZyme.count %>%
  select(colnames(orthogroups.count)[colnames(orthogroups.count) %in% metadata$strain[metadata$name == "G. tritici"]]) %>%
  filter(rowSums(across())>0)

#Filter orthogroup presence-absence matrix for orthogroups with at least one predicted CSEP/CAZyme
CSEP.orthogroups.gt <- orthogroups.count.gt[
  match(rownames(CSEP.count.gt)[which(rowSums(CSEP.count.gt) > 0)],
        rownames(orthogroups.count.gt)),]
CAZyme.orthogroups.gt <- orthogroups.count.gt[
  match(rownames(CAZyme.count.gt)[which(rowSums(CAZyme.count.gt) > 0)],
        rownames(orthogroups.count.gt)),]

#Categorise Gt orthogroups as core/accessory/specific
orthogroups.stats.gt <- data.frame(orthogroups=rownames(orthogroups.count.gt),
                                   category=NA)

progress.bar <- txtProgressBar(1, length(rownames(orthogroups.count.gt)), initial=0, char="=", style=3)
for (i in 1:length(rownames(orthogroups.count.gt))) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, i)
  
  if (length(which(orthogroups.count.gt[i,] == 0)) == (length(colnames(orthogroups.count.gt)) - 1)) {
    orthogroups.stats.gt$category[orthogroups.stats.gt$orthogroup == rownames(orthogroups.count.gt)[i]] <-
      "specific"
  }
  if (length(which(orthogroups.count.gt[i,] == 0)) == 0) {
    orthogroups.stats.gt$category[orthogroups.stats.gt$orthogroup == rownames(orthogroups.count.gt)[i]] <-
      "core"
  }
  if (length(which(orthogroups.count.gt[i,] > 0)) == (length(colnames(orthogroups.count.gt)) - 1)) {
    orthogroups.stats.gt$category[orthogroups.stats.gt$orthogroup == rownames(orthogroups.count.gt)[i]] <-
      "softcore"
  }
  if (length(which(orthogroups.count.gt[i,] == 0)) < (length(colnames(orthogroups.count.gt)) - 1) &&
      length(which(orthogroups.count.gt[i,] == 0)) > 1) {
    orthogroups.stats.gt$category[orthogroups.stats.gt$orthogroup == rownames(orthogroups.count.gt)[i]] <-
      "accessory"
  }
}
close(progress.bar)

#Make dataframe summarising Gt gene content
genes.df.gt <- data.frame(strain=rep(colnames(orthogroups.count.gt), each=4), 
                          category=rep(c("specific", "accessory", "softcore", "core"),
                                       length(colnames(orthogroups.count.gt))),
                          orthogroups=NA,
                          CSEPs=NA,
                          CAZymes=NA,
                          BGCs=NA)

#For each strain...
for (i in unique(genes.df.gt$strain)) {
  
  #For each category...
  for (j in unique(genes.df.gt$category)) {
    
    #Get number of orthogroups
    genes.df.gt$orthogroups[intersect(which(genes.df.gt$strain == i), which(genes.df.gt$category == j))] <-
      table(orthogroups.stats.gt$category[
        match(rownames(orthogroups.count.gt[orthogroups.count.gt[, i] > 0,]),
              orthogroups.stats.gt$orthogroup)])[j]
    
    #Get number of CSEPs
    genes.df.gt$CSEPs[intersect(which(genes.df.gt$strain == i), which(genes.df.gt$category == j))] <-
      table(orthogroups.stats.gt$category[
        match(rownames(CSEP.orthogroups.gt[CSEP.orthogroups.gt[, i] > 0,]),
              orthogroups.stats.gt$orthogroup)])[j]
    
    #Get number of CAZymes
    genes.df.gt$CAZymes[intersect(which(genes.df.gt$strain == i), which(genes.df.gt$category == j))] <-
      table(orthogroups.stats.gt$category[
        match(rownames(CAZyme.orthogroups.gt[CAZyme.orthogroups.gt[, i] > 0,]),
              orthogroups.stats.gt$orthogroup)])[j]
    
    #Get number of BGCs
    genes.df.gt$BGCs[intersect(which(genes.df.gt$strain == i), which(genes.df.gt$category == j))] <-
      table(bgc.num.df.gt$category[bgc.num.df.gt$strain == i])[j]
    
  }
}

#Melt for plotting
genes.df.gt <- genes.df.gt %>%
  gather("variable", "value", -strain, -category)

#Summarise non-Gt gene content
genes.df.other <- bind_rows(
  orthogroups.count %>%
    mutate(orthogroup=rownames(orthogroups.count),
           biotype=orthogroups.stats$biotype[match(orthogroup, orthogroups.stats$orthogroup)]) %>%
    filter(is.na(biotype)) %>%
    select(c("orthogroup",
             colnames(orthogroups.count)[
               colnames(orthogroups.count) %in%
                 metadata$strain[metadata$name != "G. tritici"]
             ])) %>%
    gather(strain, num, -orthogroup) %>%
    filter(num > 0) %>%
    group_by(strain) %>%
    summarise(value=n()) %>%
    mutate(variable="orthogroups"),
  
  CSEP.count %>%
    select(colnames(CSEP.count)[colnames(CSEP.count) %in% metadata$strain[metadata$name != "G. tritici"]]) %>%
    mutate(orthogroup=rownames(CSEP.count)) %>%
    gather(strain, num, -orthogroup) %>%
    filter(num > 0) %>%
    group_by(strain) %>%
    summarise(value=n()) %>%
    mutate(variable="CSEPs"),
  
  CAZyme.count %>%
    select(colnames(CAZyme.count)[colnames(CAZyme.count) %in% metadata$strain[metadata$name != "G. tritici"]]) %>%
    mutate(orthogroup=rownames(CAZyme.count)) %>%
    gather(strain, num, -orthogroup) %>%
    filter(num > 0) %>%
    group_by(strain) %>%
    summarise(value=n()) %>%
    mutate(variable="CAZymes"),
  
  bgc.num.df %>%
    filter(strain %in% metadata$strain[metadata$name != "G. tritici"]) %>%
    group_by(strain) %>%
    summarise(value=n()) %>%
    mutate(variable="BGCs")
) %>%
  mutate(category=NA)

#Combine Gt and non-Gt gene content
genes.df <- bind_rows(genes.df.gt, genes.df.other) %>%
  mutate(category=factor(category, levels=c("specific", "accessory", "softcore", "core")),
         tip=metadata$tip[match(strain, metadata$strain)])

#Add empty outgroup rows
genes.df <- bind_rows(genes.df,
                      data.frame(tip=metadata$tip[metadata$clade == "outgroup"],
                                 value=NA,
                                 variable="orthogroups"))  %>%
  mutate(variable=factor(variable, levels=c("orthogroups", "CSEPs", "CAZymes", "BGCs")))

#Make function for dynamic axis labels
addUnits <- function(n) {
  labels <- ifelse(n < 1000, n,
                   ifelse(n < 1e6, paste0(round(n/1e3), 'k')
                   )
  )
  return(labels)
}

#Plot bargraphs of number of genes
gg.gene.numbers <- ggplot(genes.df, aes(y=tip, x=value, fill=category)) +
  facet_wrap(. ~ variable, scales="free_x", nrow=1,
             labeller=labeller(variable=c(CSEPs="CSEPs", CAZymes="CAZymes", BGCs="BGCs", orthogroups="All genes"))) +
  geom_bar(stat="identity", linewidth=0.5, width=0.6) +
  geom_label(data=genes.df %>% group_by(tip, variable) %>% summarise(total=sum(na.omit(value))) %>% filter(total != 0),
             aes(x=total, y=tip, label=comma(total)), 
             position=position_stack(vjust=0.5),
             size=1.5,
             label.padding=unit(0.1, "lines"),
             label.size=0,
             alpha=0.5,
             inherit.aes=FALSE) +
  labs(x="Total numbers") +
  guides(fill=guide_legend(nrow=1)) +
  scale_x_continuous(expand=c(0, 0),
                     position="top",
                     labels=addUnits) +
  scale_fill_manual(values=c("#F5D896", "#edb945", "#e69f00", "#a58337"),
                    #c("grey93", "lightgrey", "darkgrey", "dimgrey"),
                    breaks=c("core", "softcore", "accessory", "specific"),
                    labels=c("Core", "Soft-core", "Accessory", "Specific")) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(strip.placement="outside",
        strip.text=element_text(size=6, margin=margin(0, 0, 1, 0)),
        axis.title.x=element_text(size=7, face="bold"),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position=c(-0.7, 1.04),
        legend.text=element_text(size=6, margin=margin(0, 5, 0, 0)),
        legend.key.size=unit(0.2, "cm"),
        legend.spacing.x=unit(0.1, "cm"),
        legend.title=element_blank())


## Gene copy number ##

#Remove TEs from orthogroup abundance matrix
copynum.df <- data.frame(strain=colnames(orthogroups.count),
                         as.data.frame(t(orthogroups.count))) %>%
  gather("hog", "num", -strain) %>%
  mutate(biotype=orthogroups.stats$biotype[match(hog, orthogroups.stats$orthogroup)]) %>%
  filter(is.na(biotype), num != 0) %>%
  mutate(category=factor(
    case_when(
      strain %in% metadata$strain[metadata$name == "G. tritici"] ~
        orthogroups.stats.gt$category[
          match(hog,
                orthogroups.stats.gt$orthogroup)]),
    levels=c("specific", "accessory", "softcore", "core")),
    tip=metadata$tip[match(strain, metadata$strain)])

#Plot jittered scatterplot of gene copynumber
gg.copynum <- ggplot(copynum.df, aes(x=num, y=tip, colour=category)) +
  geom_point(position="jitter", size=0.3) +
  coord_cartesian(clip="off") +
  scale_x_continuous(limits=c(1, NA),
                     breaks=pretty_breaks(),
                     expand=c(0, 0),
                     position="top") +
  scale_colour_manual(values=c("#F5D896", "#edb945", "#e69f00", "#a58337"),
                      breaks=c("core", "softcore", "accessory", "specific"),
                      labels=c("Core", "Soft-core", "Accessory", "Specific")) +
  labs(x="Gene copy\nnumber") +
  theme_minimal() +
  theme(legend.position="none",
        axis.title.x.top=element_text(size=7, vjust=-6, face="bold"),
        axis.text.x.top=element_text(size=6, vjust=-6),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid.major.x=element_line(colour="white"),
        panel.grid.minor.x=element_line(colour="white"),
        panel.grid.major.y=element_blank())


## Euler plot of high copy-number (HCN) genes ## 

#Filter for HCN genes
hcn.df <- copynum.df %>%
  filter(hog %in% unique(copynum.df %>%
                           filter(num > 10) %>%
                           pull(hog))) %>%
  filter(num != 0) %>%
  mutate(clade=metadata$clade[match(strain, metadata$strain)])

#Make list for collecting overlapping groups
hcn.groups <- c(
  "Gh"=0,
  "Ga"=0,
  "GtB"=0,
  "GtA"=0,
  
  "Ga&Gh"=0,
  "Gh&GtB"=0,
  "Gh&GtA"=0,
  "Ga&GtB"=0,
  "Ga&GtA"=0,
  "GtA&GtB"=0,
  
  "Ga&Gh&GtB"=0,
  "Ga&Gh&GtA"=0,
  "Gh&GtA&GtB"=0,
  "Ga&GtA&GtB"=0,
  
  "Ga&Gh&GtA&GtB"=0
)

#For each orthogroup...
for (ortho in unique(hcn.df$hog)) {
  
  #Get clades that orthogroup is present in
  clade <- hcn.df %>%
    filter(hog == ortho) %>%
    pull(clade)
  
  if (length(clade) > 1) {
    
    clade <- paste0(sort(unique(clade)), collapse="&")
    
  }
  
  #Add count to results list
  hcn.groups[clade] <- hcn.groups[clade] + 1
  
}

#Plot euler
set.seed(1)
selection.euler <- euler(hcn.groups, shape="ellipse")
eulerr_options(padding=unit(0.5, "pt"))
euler <- plot(selection.euler, 
              fills=list(fill=c(rep("white", 10), "#87AE88", rep("white", 3), "#edb945")),
              labels=list(cex=0.7,
                          col="black",
                          padding=unit(c(1, 1), "pt")),
              edges=list(col="dimgrey", 
                         #lty=c("solid", "dotted", "dotted"), 
                         lwd=0.5),
              quantities=list(cex=0.7,
                              col="black"))

#Convert to grob
euler.grob <- as.grob(euler)

#Plot dummy legend
gg.dummy <- ggplot(
  data=data.frame(
    x=0, y=0,
    colour=c("Present in GtA",
             "Absent in GtA")),
  aes(x=0, y=0, colour=colour)) +
  geom_point(size=3) +
  scale_colour_manual(values=c("#87AE88", "#edb945")) +
  labs(colour="Distribution of genes with\ncopy number outliers (>10)") +
  theme_minimal() +
  theme(legend.title=element_text(size=6, face="bold"),
        legend.key.size=unit(0.4, "cm"),
        legend.position=c(0.22, 0.85),
        legend.text=element_text(size=6))

euler.legend <- get_legend(gg.dummy)

#Combine plot with legend
gg.euler <- as_ggplot(euler.legend) +
  annotation_custom(euler.grob, xmin=0.25, xmax=Inf, ymin=-Inf, ymax=0.95)

#Write table of HCNs for GO enrichment analysis
write.table(file="copynum_hogs.tsv", data.frame(unique(hcn.df$hog)),
            col.names=FALSE, quote=FALSE, row.names=FALSE)


## Gt pangenome accumulation curves ##

#Function from https://github.com/SioStef/panplots
panplots <- function(data, curve = "pan", iterations = 100) {
  
  nr_rows <- nrow(data);
  nr_iterations <- iterations; #the number of iterations (100 by default) 
  #create empty matrix to store temp results
  temp <- matrix(data=NA,nrow=nr_rows,ncol=nr_iterations)
  
  if(curve == "core") {
    ## compute core_genome_accumulation_curve data for the number of iterations
    for(times in 1:nr_iterations){
      # random sampling of genomes
      for (i in 1: nr_rows){
        t=data[sample(nr_rows, i), ,drop=F]
        temp[i,times]=length(which(colSums(t) == i))
      }
    }
  } 
  
  if(curve == "pan") {
    ## compute gene_cluster_accumulation_curve data for the number of iterations
    for(times in 1:nr_iterations){
      # random sampling of genomes
      for (i in 1: nr_rows){
        t=data[sample(nr_rows, i), ,drop=F]
        temp[i,times]=length(which(colSums(t) > 0))
      }
    }
  } 
  
  if(curve == "uniq") {
    ## compute unique gene_cluster_accumulation_curve data for the number of iterations
    for(times in 1:nr_iterations){
      # random sampling of genomes
      for (i in 2: nr_rows){
        t=data[sample(nr_rows, i), ,drop=F]
        temp[i,times]=length(which(colSums(t) == 1))
      }
    }
  }
  
  # summerize permutation results using "matrixStats" library
  summary <- data.frame(genomes=c(1:nr_rows)) 
  summary$mean=rowMeans2(temp[,c(-1)])
  summary$sd=rowSds(temp[,c(-1)])
  summary$group=deparse(substitute(data))
  return(summary)
}

#Format abundance matrix
gt.matrix <- orthogroups.count.gt
gt.matrix[gt.matrix>0] <- 1
gt.matrix <- t(gt.matrix)

#Get data for accumulation curves
accumulation.df <- rbind(data.frame(panplots(gt.matrix, curve="core"), type="core"), 
                         data.frame(panplots(gt.matrix, curve="pan"), type="pan"))

core.num <- accumulation.df$mean[
  intersect(which(accumulation.df$type == "core"), which(accumulation.df$genomes == 5))
]

pan.num <- accumulation.df$mean[
  intersect(which(accumulation.df$type == "pan"), which(accumulation.df$genomes == 5))
]

#Plot accumulation curves
gg.accumulation <- ggplot(accumulation.df, aes(x=genomes, y=mean)) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, group=type),
              alpha=0.2, show.legend=FALSE) +
  geom_line(aes(linetype=type), linewidth=0.5) +
  # geom_smooth(method="lm", formula = y ~ I(log(x)),
  #             colour="black", se=FALSE,
  #             linewidth=0.5) +
  geom_point(size=0.5) +
  annotate(geom="text",
           x=4.8,
           y=core.num+600,
           label=paste0("Core genes:\n", prettyNum(core.num, big.mark=",")),
           fontface="bold",
           size=2) +
  annotate(geom="text",
           x=4.8,
           y=pan.num-600,
           label=paste0("Pan genes:\n", prettyNum(pan.num, big.mark=",")),
           fontface="bold",
           size=2) +
  geom_segment(x=4.9, xend=5, y=core.num+300, yend=core.num+100,
               linewidth=0.3,
               arrow=arrow(length=unit(0.02, "npc"))) +
  geom_segment(x=4.9, xend=5, y=pan.num-300, yend=pan.num-100,
               linewidth=0.3,
               arrow=arrow(length=unit(0.02, "npc"))) +
  labs(x="Number of Gt genomes sampled", y="Number of genes") +
  scale_y_continuous(labels=addUnits) +
  scale_linetype(breaks=c("pan", "core"),
                 labels=c("Pangenome", "Core genome")) +
  theme_minimal() +
  theme(legend.position=c(0.15, 0.9),
        legend.title=element_blank(),
        panel.grid.minor.x=element_blank(),
        legend.key.size=unit(8, "pt"),
        legend.text=element_text(size=7, margin=margin(0, 5, 0, 0)),
        axis.title=element_text(size=7),
        axis.text=element_text(size=5))


## COMBINE PLOTS ##

#Add labels to species tree
gg.tree2 <- gg.tree +
  geom_highlight(node=MRCA(tree, metadata$tip[metadata$name == "G. tritici" & metadata$own == "Y"]),
                 fill=NA,
                 colour="black",
                 linetype="dashed",
                 alpha=0, extend=0.105,
                 show.legend=FALSE) +
  geom_rect(xmin=0.001, xmax=0.043, ymin=6.25, ymax=11.3,
            fill=NA,
            colour="black",
            linewidth=0.3) +
  annotate(geom="text",
           label="Gene variance (PERMANOVA)",
           x=0.022,
           y=10.5,
           vjust=-1,
           fontface="bold",
           size=2.5) +
  annotate(geom="text",
           label="All genes",
           x=0.01,
           y=10,
           fontface="bold",
           size=2) +
  annotate(geom="label",
           label=paste0("Phylogeny: ", round(sum(permanova.orthogroup$R2[1:2]) * 100), "%\n",
                        "Lifestyle: ", round(sum(permanova$R2[3]) * 100), "%"),
           x=0.03,
           y=10,
           size=2) +
  annotate(geom="text",
           label="CSEPs",
           x=0.01,
           y=9,
           fontface="bold",
           size=2) +
  annotate(geom="label",
           label=paste0("Phylogeny: ", round(sum(permanova.CSEP$R2[1:2]) * 100), "%\n",
                        "Lifestyle: ", round(sum(permanova.CSEP$R2[3]) * 100), "%"),
           x=0.03,
           y=9,
           size=2) +
  annotate(geom="text",
           label="CAZymes",
           x=0.01,
           y=8,
           fontface="bold",
           size=2) +
  annotate(geom="label",
           label=paste0("Phylogeny: ", round(sum(permanova.CAZyme$R2[1:2]) * 100), "%\n",
                        "Lifestyle: ", round(sum(permanova.CAZyme$R2[3]) * 100), "%"),
           x=0.03,
           y=8,
           size=2) +
  annotate(geom="text",
           label="BGCs",
           x=0.01,
           y=7,
           fontface="bold",
           size=2) +
  annotate(geom="label",
           label=paste0("Phylogeny: ", round(sum(permanova.BGC$R2[1:2]) * 100), "%\n",
                        "Lifestyle: ", round(sum(permanova.BGC$R2[3]) * 100), "%"),
           x=0.03,
           y=7,
           size=2)

#Set separate legends for composite plots
options("aplot_guides"="keep")

#Combine tree with gene content summaries
gene.content.grob <- as.grob(
  gg.gene.numbers %>%
    insert_right(gg.copynum, width=0.4) %>%
    insert_left(gg.tree2, width=2.7)
)

#Write to file
pdf(file=paste0("R:/GaeumannomycesGenomics/07_comparative_genomics/gaeumannomyces_genes-",
                Sys.Date(), ".pdf"),
    height=6, width=7)
ggarrange(gene.content.grob,
          ggarrange(gg.accumulation, gg.euler, labels=c("", "c")),
          nrow=2, heights=c(2, 1), labels="auto")
dev.off()


## Abundance matrix plots ##

#Remove other Gt tip from tree
skeleton.tree <- drop.tip(tree, c(outgroups, "GCF_000145635.1_Gae_graminis_V2_protein.faa_XP"))

#Plot bare tree
gg.skeleton.tree <- ggtree(skeleton.tree, branch.length="none") %<+% metadata +
  xlim(0, 30) +
  geom_tiplab(aes(label=new.label),
              align=TRUE,
              offset=0.5,
              size=2,
              parse=TRUE) +
  theme(plot.title=element_text(face="bold"))


# CSEPs #

#Read in PHI-base data
phibase.df <- read.csv("R:/GaeumannomycesGenomics/03_functional_annotation/phi-base_current.csv")

#Summarise number of CSEPs with PHI-base names 
CSEP.abundance.df <- CSEP.count %>%
  filter(rowSums(across())>0) %>%
  rownames_to_column(var="orthogroup") %>%
  gather(strain, num, -orthogroup) %>%
  mutate(PHI.base_entry=sub("_.*", "",
                            orthogroups.stats$PHI.base_entry[match(orthogroup, orthogroups.stats$orthogroup)]),
         CSEP_name=orthogroups.stats$CSEP_name[match(orthogroup, orthogroups.stats$orthogroup)],
         name_label=paste0(sub("N0\\.", "", orthogroup), "\n", CSEP_name)) %>%
  filter(num > 0, !is.na(PHI.base_entry)) %>%
  mutate(phenotype=NA)

#Add mutant phenotype groups to dataframe
for (i in 1:length(unique(CSEP.abundance.df$PHI.base_entry))) {
  
  id <- unlist(strsplit(unique(CSEP.abundance.df$PHI.base_entry)[i], ","))
  
  CSEP.abundance.df$phenotype[CSEP.abundance.df$PHI.base_entry == unique(CSEP.abundance.df$PHI.base_entry)[i]] <- 
    paste(unique(phibase.df$Phenotype_of_mutant[which(phibase.df$PHI_MolConn_ID %in% id)]), collapse=",")
  
}

#Split rows with more than one substrate and add tree tip labels
CSEP.abundance.df <- CSEP.abundance.df %>%
  separate_rows(sep=",", phenotype) %>%
  group_by(phenotype) %>%
  complete(strain, name_label) %>%
  mutate(tiplab=metadata$tip[match(strain, metadata$strain)],
         phenotype=sub("Effector \\(plant avirulence determinant\\)", "Effector", 
                       sub("Loss of pathogenicity", "Loss\npath", str_to_sentence(phenotype)))) %>%
  select(tiplab, everything())

#Plot grid of CSEPs
gg.cseps <- ggplot(CSEP.abundance.df,
                   aes(x=name_label, y=tiplab, fill=num)) +
  facet_grid(. ~ phenotype, scales="free", space="free") +
  geom_tile(color="grey", linewidth=0.1) +
  scale_fill_gradient(low="#F5D896", high="#a58337",
                      breaks=pretty_breaks(),
                      na.value="white",
                      guide=guide_colourbar(
                        title="CSEP copy-number",
                        title.position="left",
                        direction="horizontal",
                        title.hjust=0,
                        title.vjust=0.8)) +
  theme_minimal() +
  theme(legend.position="top",
        legend.direction="horizontal",
        legend.margin=margin(5, 0, -10, 0),
        legend.title=element_text(size=5, face="bold",
                                  margin=margin(0, 10, 0, 0)),
        legend.text=element_text(size=3, margin=margin(0, 3, 0, 0)),
        legend.key.size=unit(6, "pt"),
        strip.text=element_text(size=4, face="bold",
                                margin=margin(2, 0, 2, 0)),
        panel.spacing=unit(0.15, "lines"),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5,
                                 size=3, margin=margin(-2, 0, 0, 0)))


# CAZymes #

#Summarise number of CAZymes with known PCWDE substrates
CAZyme.abundance.df <- CAZyme.count %>%
  filter(rowSums(across())>0) %>%
  rownames_to_column(var="orthogroup") %>%
  gather(strain, num, -orthogroup) %>%
  mutate(CAZyme_family=orthogroups.stats$CAZyme_family[match(orthogroup, orthogroups.stats$orthogroup)],
         CAZyme_family=gsub("\\+[[:digit:]].*$", "\\1\\",
                            gsub("\\+[[:digit:]].*\\+", "\\1\\+",
                                 sub(",.*", "",
                                     sub("^-,", "", CAZyme_family))))) %>%
  filter(num > 0) %>%
  group_by(strain, CAZyme_family) %>%
  summarise(count=n()) %>%
  mutate(substrate=NA)

#Get substrate groupings
substrates.df <- read.csv("R:/GaeumannomycesGenomics/07_comparative_genomics/cazyme_substrates.csv")

#Add substrates to dataframe
for (i in 1:length(unique(substrates.df$CAZy.Family))) {
  
  substrates <- substrates.df$Substrate[substrates.df$CAZy.Family == unique(substrates.df$CAZy.Family)[i]]
  
  if (length(substrates) > 0) {
    
    CAZyme.abundance.df$substrate[grep(paste0("\\b", unique(substrates.df$CAZy.Family)[i], "\\b"),
                                       CAZyme.abundance.df$CAZyme_family)] <- paste(substrates, collapse=",") 
    
  }
  
}

#Split rows with more than one substrate and add tree tip labels
CAZyme.abundance.df<- CAZyme.abundance.df %>%
  separate_rows(sep=",", substrate) %>%
  filter(!is.na(substrate)) %>%
  group_by(substrate) %>%
  complete(strain, CAZyme_family) %>%
  mutate(tiplab=metadata$tip[match(strain, metadata$strain)],
         substrate=sub("Cutin", "Cu", substrate)) %>%
  select(tiplab, everything())

#Plot grid of CAZymes
gg.cazymes <- ggplot(CAZyme.abundance.df,
                     aes(x=CAZyme_family, y=tiplab, fill=count)) +
  facet_grid(. ~ substrate, scales="free", space="free") +
  geom_tile(color="grey", linewidth=0.1) +
  # scale_x_discrete(labels=bold_labels_cseps) +
  scale_fill_gradient(low="#F5D896", high="#a58337",
                      breaks=pretty_breaks(),
                      na.value="white",
                      guide=guide_colourbar(
                        title="Number of CAZymes",
                        title.position="left",
                        direction="horizontal",
                        title.hjust=0,
                        title.vjust=0.8)) +
  theme_minimal() +
  theme(legend.position="top",
        legend.direction="horizontal",
        legend.margin=margin(5, 0, -10, 0),
        legend.title=element_text(size=5, face="bold",
                                  margin=margin(0, 10, 0, 0)),
        legend.text=element_text(size=3, margin=margin(0, 3, 0, 0)),
        legend.key.size=unit(6, "pt"),
        strip.text=element_text(size=4, face="bold",
                                margin=margin(2, 0, 2, 0)),
        panel.spacing=unit(0.15, "lines"),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5,
                                 size=3, margin=margin(-2, 0, 0, 0)))

# BGCS #

#Summarise number of BGCs for each class
bgc.classes.df <- bgc.clusters.df %>%
  rename(BGC="X.BGC.Name") %>%
  inner_join(bgc.df)

bgc.abundance.df <- data.frame(t(bgc.abundance.mat %>%
                                   column_to_rownames(var="strain"))) %>%
  rownames_to_column(var="cluster") %>%
  gather(strain, num, -cluster) %>%
  mutate(num=as.numeric(sub(0, NA, num)),
         class=bgc.classes.df$BiG.SCAPE.class[match(cluster, bgc.classes.df$Family.Number)],
         tiplab=metadata$tip[match(sub("\\.", "-", strain), metadata$strain)]) %>%
  select(tiplab, everything())

#Make function for plot labels
bgc_label <- function(breaks) {
  labels <- paste0("BGC", breaks)
}

#Plot grid of BGCs
gg.bgcs <- ggplot(bgc.abundance.df,
                  aes(x=cluster, y=tiplab, fill=num)) +
  facet_grid(. ~ class, scales="free", space="free") +
  geom_tile(color="grey", linewidth=0.1) +
  scale_fill_gradient(low="#F5D896", high="#a58337",
                      breaks=pretty_breaks(),
                      na.value="white",
                      guide=guide_colourbar(
                        title="BGC copy-number",
                        title.position="left",
                        direction="horizontal",
                        title.hjust=0,
                        title.vjust=0.8)) +
  scale_x_discrete(labels=bgc_label) +
  theme_minimal() +
  theme(legend.position="top",
        legend.direction="horizontal",
        legend.margin=margin(5, 0, -10, 0),
        legend.title=element_text(size=5, face="bold",
                                  margin=margin(0, 10, 0, 0)),
        legend.text=element_text(size=3, margin=margin(0, 3, 0, 0)),
        legend.key.size=unit(6, "pt"),
        strip.text=element_text(size=4, face="bold",
                                margin=margin(2, 0, 2, 0)),
        panel.spacing=unit(0.15, "lines"),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5,
                                 size=3, margin=margin(-2, 0, 0, 0)))


#Write to file
pdf(file=paste0("R:/GaeumannomycesGenomics/07_comparative_genomics/abundances-1-", Sys.Date(), ".pdf"),
    height=4.2, width=7)
ggarrange(as.grob(gg.cseps %>% insert_left(gg.skeleton.tree, width=0.27)),
          as.grob(gg.cazymes %>% insert_left(gg.skeleton.tree, width=0.27)),
          labels="auto", nrow=2, heights=c(1.1, 1))
dev.off()

pdf(file=paste0("R:/GaeumannomycesGenomics/07_comparative_genomics/abundances-2-", Sys.Date(), ".pdf"),
    height=1.8, width=7)
ggarrange(as.grob(gg.bgcs %>% insert_left(gg.skeleton.tree, width=0.27)),
          labels="c", nrow=1)
dev.off()


## Repeat classification ##

#Read in eirepeat results
repeats.dir <- "S:/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Analysis/"

for (strain in metadata$strain[metadata$own == "Y"]) {
  
  repeats <- rbind(
    read.table(paste0(repeats.dir, strain, "/Analysis/eirepeat-1.1.0/output/RepeatMasker_interspersed_repeatmodeler/genome.fa.out"),
               skip=3, header=FALSE, fill=TRUE,
               col.names=c("sw.score", "perc.div", "perc.del", "perc.ins", "sequence",
                           "sequence.begin", "sequence.end", "sequence.left", "C",
                           "matching.repeat", "repeatclass.family", "repeat.begin",
                           "repeat.end", "repeat.left", "id", "asterisk")),
    read.table(paste0(repeats.dir, strain, "/Analysis/eirepeat-1.1.0/output/RepeatMasker_interspersed/genome.fa.out"),
               skip=3, header=FALSE, fill=TRUE,
               col.names=c("sw.score", "perc.div", "perc.del", "perc.ins", "sequence",
                           "sequence.begin", "sequence.end", "sequence.left", "C",
                           "matching.repeat", "repeatclass.family", "repeat.begin",
                           "repeat.end", "repeat.left", "id", "asterisk")))
  
  repeats.index <- repeats %>%
    select(matching.repeat, repeatclass.family) %>%
    distinct()
  
  #Summarise number of repeats for each class
  repeats.gff <- read.table(paste0(repeats.dir, strain, "/Analysis/eirepeat-1.1.0/output/all_interspersed_repeats.gff3")) %>%
    filter(V3 == "match") %>%
    mutate(motif=sub("^.*Motif:", "", V9)) %>%
    mutate(repeatclass.family=repeats.index$repeatclass.family[match(motif, repeats.index$matching.repeat)]) %>%
    separate(repeatclass.family, c("class", "family"), sep="/") %>%
    group_by(class, family) %>%
    summarise(num=n()) %>%
    mutate(tip=metadata$tip[metadata$strain == strain])
  
  assign(paste0(strain, ".repeats"), repeats.gff)
  
}

#Combine repeat annotations for all strains
repeats.df <- bind_rows(mget(ls(pattern="\\.repeats"))) %>%
  filter(class != "ARTEFACT")
#Add empty row for outgroup
repeats.df <- bind_rows(repeats.df, data.frame(tip=metadata$tip[metadata$clade == "outgroup"]))

#Plot bargraph of number and classifcation of repeats
gg.repeats <- ggplot(repeats.df, aes(y=tip, x=num, fill=class)) +
  geom_bar(stat="identity", linewidth=0.5, width=0.6) +
  labs(x="Number of TEs") +
  scale_x_continuous(expand=c(0, 0),
                     position="top",
                     labels=addUnits) +
  scale_fill_manual(values=c("#88CCEE", "#44AA99", "#117733", "#332288",
                                      "#DDCC77", "#999933","#CC6677",
                                      "#882255", "#AA4499", "dimgrey"),
                                      labels=sub("_", " ", sort(unique(repeats.df$class)))) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(strip.placement="outside",
        strip.text=element_text(size=6, margin=margin(0, 0, 1, 0)),
        axis.title.x=element_text(size=8),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x.top=element_text(size=7),
        panel.grid.major.x=element_line(colour="white"),
        panel.grid.minor.x=element_line(colour="white"),
        panel.grid.major.y=element_blank(),
        legend.position=c(-2, 0.8),
        legend.text=element_text(size=8, margin=margin(0, 5, 0, 0)),
        legend.key.size=unit(0.2, "cm"),
        legend.spacing.x=unit(0.1, "cm"),
        legend.title=element_blank())

pdf(file=paste0("R:/GaeumannomycesGenomics/07_comparative_genomics/gaeumannomyces_TEs-",
                Sys.Date(), ".pdf"),
    height=4, width=7)
gg.repeats %>% 
  insert_left(gg.tree, width=2.7)
dev.off()
