#Written in R v4.3.1
library(tidyverse)    #2.0.0
library(ape)          #5.7-1
library(cowplot)      #1.1.1
library(eulerr)       #7.0.0
library(ggforce)      #0.4.1
library(ggh4x)        #0.2.6
library(gggenomes)    #0.9.12.9000
library(ggmsa)        #1.6.0
library(ggnewscale)   #0.4.9
library(ggplotify)    #0.1.2
library(ggpubr)       #0.6.0
library(ggrepel)      #0.9.3
library(ggtree)       #3.9.1
library(multcompView) #0.1-9
library(regioneR)     #1.32.0
library(rstatix)      #0.7.2
library(rtracklayer)  #1.60.1
library(scales)       #1.2.1
library(seqmagick)    #0.1.6


#Directory paths
dir.genephylo <- "R:/GaeumannomycesGenomics/04_phylogenetic_classification/"
dir.phylo <- "R:/GaeumannomycesGenomics/05_phylogenomics/"
dir.synteny <- "R:/GaeumannomycesGenomics/06_synteny/"
dir.comp <- "R:/GaeumannomycesGenomics/07_comparative_genomics/"


#Create functions to fix ordering of clades and strains in plots
plot_orders <- function(df) {
  df %>%
    mutate(clade=factor(metadata$clade[match(new.strain, metadata$new.strain)],
                        levels=c("GtB", "GtA", "Ga", "Gh")),
           clade.label=recode_factor(clade,
                                     "GtB"='paste(italic("Gt"), "B")',
                                     "GtA"='paste(italic("Gt"), "A")',
                                     "Ga"='paste(italic("Ga"))',
                                     "Gh"='paste(italic("Gh"))'),
           new.strain=factor(new.strain,
                             levels=c("Gt-LH10", "Gt-4e", "Gt-23d", "Gt-8d", "Gt-19d1",
                                      "Ga-3aA1", "Ga-CB1", "Gh-2C17", "Gh-1B17")),
           new.strain.label=recode_factor(new.strain,
                                          "Gt-LH10"='paste("Gt-LH10")',
                                          "Gt-4e"='paste("Gt-4e")',
                                          "Gt-23d"='paste("Gt-23d")',
                                          "Gt-8d"='paste("Gt-8d")',
                                          "Gt-19d1"='paste("Gt-19d1")',
                                          "Ga-3aA1"='paste("Ga-3aA1")',
                                          "Ga-CB1"='paste("Ga-CB1")',
                                          "Gh-2C17"='paste("Gh-2C17")',
                                          "Gh-1B17"='paste("Gh-1B17")'))
}

plot_orders_rev <- function(df) {
  df %>%
    mutate(clade=factor(metadata$clade[match(new.strain, metadata$new.strain)],
                        levels=c("Gh", "Ga", "GtA", "GtB")),
           clade.label=recode_factor(clade,
                                     "Gh"='paste(italic("Gh"))',
                                     "Ga"='paste(italic("Ga"))',
                                     "GtA"='paste(italic("Gt"), "A")',
                                     "GtB"='paste(italic("Gt"), "B")'),
           new.strain=factor(new.strain,
                             levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                                      "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10")),
           new.strain.label=recode_factor(new.strain,
                                          "Gh-1B17"='paste("Gh-1B17")',
                                          "Gh-2C17"='paste("Gh-2C17")',
                                          "Ga-CB1"='paste("Ga-CB1")',
                                          "Ga-3aA1"='paste("Ga-3aA1")',
                                          "Gt-19d1"='paste("Gt-19d1")',
                                          "Gt-8d"='paste("Gt-8d")',
                                          "Gt-23d"='paste("Gt-23d")',
                                          "Gt-4e"='paste("Gt-4e")',
                                          "Gt-LH10"='paste("Gt-LH10")'))
}


#Make list with strains and filenames
strains <- list(strain=c("Gt-19d1", "Gt-8d", "Gt-23d", "Gt14LH10", "Gt-4e",
                         "Gt-CB1", "Gt-3aA1", "Gh-2C17", "Gh-1B17"),
                file1=c("Gt17LH19d1", "Gt17LH48d", "Gt23d", "Gt14LH10", "Gt4e",
                        "GtCB1", "PG3aA1", "NZ1292C17", "Gh1B17"),
                file2=c("Gt-19d1", "Gt-8d", "Gt-23d", "Gt14LH10", "Gt-4e",
                        "Gt-CB1", "Gt-3aA1", "Gh-2C17", "Gh-1B17"))

#Read in synteny information inferred from GENESPACE
synteny <- read.csv(paste0(dir.synteny, "pseudochromosomes.tsv"),
                    sep="\t", header=FALSE) %>%
  dplyr::rename(strain="V1", full.contig="V2",
                pseudochromosome="V3", inversion="V4") %>%
  mutate(contig=sub(".*_", "", full.contig),
         pseudochromosome.abb=sub("-.*", "", pseudochromosome),
         pseudochromosome.col=ifelse(grepl("B", pseudochromosome.abb),
                                     "mixed", pseudochromosome.abb))

#Read in metadata
metadata <- read.csv(paste0(dir.phylo, "raxmlng/metadata.csv"))

#Load orthogroup data
load(paste0(dir.comp, "orthogroup-matrices-2023-07-25.RData"))

#Remove TEs from orthogroup matrix
orthogroups <- orthogroups %>%
  rownames_to_column(var="orthogroup") %>%
  mutate(biotype=orthogroups.stats$biotype[match(orthogroup, orthogroups.stats$orthogroup)]) %>%
  filter(is.na(biotype)) %>%
  column_to_rownames(var="orthogroup") %>%
  select(-biotype)


################################################################################
###################### IDENTIFY DUPLICATE GENE LOCATIONS #######################
################################################################################

#Make empty list for duplicate location results
duplicates.stats <- list()

#For each strain...
for (strain in colnames(orthogroups)) {
  
  #Print progress
  message(strain)
  
  #Read in GENESPACE pangenes
  pangenes <- read.csv(paste0("S:/genespace_gaeumannomyces_230510/pangenes/",
                              sub("-", "", strain), "_pangenes.txt"), sep="\t") %>%
    filter(genome == sub("-", "", strain)) %>%
    mutate(contig=sub(".*_", "", chr)) %>%
    mutate(across(chr, ~ with(synteny, pseudochromosome[match(.x, sub("-", "", full.contig))]))) %>%
    mutate(chr2=sub("-.*$", "", chr)) %>%
    filter(!is.na(chr2))
  
  #Empty lists for different duplicate categories
  tandem.results <- list()
  intra.results <- list()
  inter.results <- list()
  mixed.results <- list()
  
  #For each orthogroup..
  progress.bar <- txtProgressBar(1, nrow(orthogroups), initial=0, char="=", style=3)
  for (i in 1:nrow(orthogroups)) {
    
    #Update progress bar
    setTxtProgressBar(progress.bar, i)
    
    #Get genes
    genes <- unlist(strsplit(orthogroups[i, strain], ", "))
    
    #If copy-number is greater than 1...
    if (length(genes) > 1) {
      
      #Get duplicates and arrange by position
      duplicates <- pangenes %>%
        filter(id %in% genes) %>%
        arrange(ord) %>%
        mutate(strain=strain)
      
      #Get tandems (genes next to eachother)
      tandems <- duplicates[c(which(diff(duplicates$ord) == 1),
                              which(diff(duplicates$ord) == 1) + 1),]
      
      if (nrow(tandems) > 0) {
        
        duplicates.stats[[length(duplicates.stats) + 1]] <- 
          data.frame(duplicate="tandem",
                     num=nrow(tandems), 
                     hog=rownames(orthogroups)[i],
                     strain=strain)
        
        tandem.results[[rownames(orthogroups)[i]]] <- tandems
        
      }
      
      #Summarise number of duplicates by pseudochromsome
      duplicates.num <- duplicates %>%
        group_by(chr2) %>%
        summarise(num=n()) %>%
        pull(num)
      
      #Count multiple duplicates on more than one pseudochromosome as mixed
      if (length(duplicates.num) > 1 && any(duplicates.num != 1)) {
        
        duplicates.stats[[length(duplicates.stats) + 1]] <-  
          data.frame(duplicate="mixed",
                     num=sum(duplicates.num), 
                     hog=rownames(orthogroups)[i],
                     strain=strain)
        
        mixed.results[[rownames(orthogroups)[i]]] <- duplicates
        
      }
      
      #Count single duplicates on more than one pseudochromosome as inter-chromosomal
      if (length(duplicates.num) > 1 && all(duplicates.num == 1)) {
        
        duplicates.stats[[length(duplicates.stats) + 1]] <-  
          data.frame(duplicate="inter-chromosomal",
                     num=sum(duplicates.num), 
                     hog=rownames(orthogroups)[i],
                     strain=strain)
        
        inter.results[[rownames(orthogroups)[i]]] <- duplicates
        
      }
      
      #Count duplicates within one pseudochromosome as intra-chromosomal
      if (length(duplicates.num) == 1) {
        
        duplicates.stats[[length(duplicates.stats) + 1]] <- 
          data.frame(duplicate="intra-chromosomal",
                     num=duplicates.num, 
                     hog=rownames(orthogroups)[i],
                     strain=strain)
        
        intra.results[[rownames(orthogroups)[i]]] <- duplicates
        
      }
      
    }
    
  } 
  close(progress.bar)
  
  df.tandems <- bind_rows(tandem.results) %>%
    distinct()
  
  assign(paste0("df.tandems.", strains$strain[strains$file2 == strain]), df.tandems)
  
}

#Combine duplicate results for all strains
df.duplicates.stats <- bind_rows(duplicates.stats) %>%
  filter(duplicate != "tandem") %>%
  mutate(duplicate=sub("mixed", "inter-chromosomal", duplicate),
         new.strain=metadata$new.strain[match(strain, metadata$strain)])

#Set plot order
df.duplicates.stats <- plot_orders_rev(df.duplicates.stats)

#Plot boxplots of duplicate set size
gg.duplicates <- ggplot(df.duplicates.stats,
                        aes(x=new.strain, y=num, colour=duplicate)) +
  facet_grid(~clade.label, scales="free_x", space="free", switch="x",
             labeller=label_parsed) +
  geom_boxplot(width=0.3,
               linewidth=0.3,
               outlier.size=0.3,
               position=position_dodge(width=0.5)) +
  labs(x=NULL, y="Duplicate gene set size") +
  scale_colour_manual(values=c("#4B854C", "#C3D6C3"),
                      labels=c("inter-chromosomal \n(including mixed)",
                               "intra-chromosomal\n(including tandems)",
                               "mixed")) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.2)),
                     limits=c(2, NA)) +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(strip.placement="outside",
        strip.text=element_text(size=7, face="bold"),
        legend.position="top",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=7, margin=margin(0, 5, 0, 0)),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=7),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=7),
        panel.grid.major.x=element_blank(),
        panel.spacing=unit(1, "lines"),
        plot.margin=margin(5.5, 5.5, 0, 5.5))


################################################################################
########################## FORMAT DATA FOR IDEOGRAMS ###########################
################################################################################

#Get predicted CSEP genes
CSEP.genes <- orthogroups[match(orthogroups.stats$orthogroup[
  !is.na(orthogroups.stats$CSEP)], rownames(orthogroups)),
]

#Make dataframe with orthogroup copy-number
copynum.df <- data.frame(strain=colnames(orthogroups.count),
                         as.data.frame(t(orthogroups.count))) %>%
  gather("hog", "num", -strain) %>%
  mutate(biotype=orthogroups.stats$biotype[
    match(hog, orthogroups.stats$orthogroup)
  ]) %>%
  filter(is.na(biotype))

#Filter for high copy-number (HCN) orthogroups
hcn.df <- copynum.df %>%
  filter(hog %in% unique(copynum.df %>%
                           filter(num > 10) %>%
                           pull(hog))) %>%
  filter(num > 0) %>%
  mutate(clade=metadata$clade[match(strain, metadata$strain)])

#Get all genes for orthogroups
all.ortho.genes <- orthogroups %>%
  rownames_to_column("hog") %>%
  mutate(biotype=orthogroups.stats$biotype[
    match(hog, orthogroups.stats$orthogroup)
  ]) %>%
  filter(is.na(biotype)) %>%
  gather(strain, gene, -hog) %>%
  separate_rows(sep=", ", gene)

#Filter for genes in HCN orthogroups
hcn.genes <- all.ortho.genes %>%
  filter(hog %in% unique(hcn.df$hog))

#Get Starship elements
starship.elements <- 
  read.csv(paste0(dir.comp, "starfish/elementFinder/gaeumannomyces.elements.bed"),
           sep="\t", header=FALSE) %>%
  dplyr::rename(seqnames="V1", start="V2", end="V3", gene="V4", category="V5",
                strain="V6", elementID="V7", flank="V8") %>%
  mutate(strain=sub("_.*", "", seqnames),
         seqnames=sub(".*_", "", seqnames))

#For each strain...
for (i in 1:length(strains$strain)) {
  
  #Print progress
  message(strains$strain[i])
  
  #Read in the header of the gff3 annotation
  header <- readLines(
    paste0("S://CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/",
           strains$file2[i], "/", strains$file2[i], "_EIv1.release.gff3"),
    n=500
  )
  
  #Pull the sequence regions
  header <- sub("   ", " ",
                header[grepl(header,
                             pattern="##sequence-region")])
  
  #Create a dataframe with the start and end positions of the sequences
  df <- data.frame(do.call(rbind, strsplit(header, split=" ")))
  df[,3] <- as.numeric(as.character(df[,3]))
  df[,4] <- as.numeric(as.character(df[,4]))
  
  #Import the whole annotation
  gff <- import(
    paste0("S://CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/",
           strains$file2[i], "/", strains$file2[i], "_EIv1.release.gff3")
  )
  
  #Add sequence lengths
  seqlengths(gff) <- 
    df$X4[match(names(seqlengths(gff)), df$X2)] -
    df$X3[match(names(seqlengths(gff)), df$X2)]
  
  #Summarise the fragments with annotated genes and format for plotting
  fragments <- GRanges(seqinfo(gff))
  names(fragments) <- sub(".*_", "", names(fragments))
  df.fragments <- as.data.frame(fragments)
  df.fragments$seqnames <- sub(".*_", "", df.fragments$seqnames)
  df.fragments$seqnames <- factor(
    df.fragments$seqnames,
    levels=df.fragments$seqnames[order(df.fragments$width, decreasing=FALSE)]
  )
  
  #Adapt start and end positions to plot in-line
  df.fragments.cumulative <- df.fragments %>%
    dplyr::slice(match(synteny$contig[synteny$strain == strains$file2[i]],
                       df.fragments$seqnames)) %>%
    mutate(colour=synteny$pseudochromosome.col[synteny$strain == strains$file2[i]],
           end2=cumsum(end),
           start2=lag(end2, default=0)+1,
           strain=strains$strain[i],
           new.strain=metadata$new.strain[match(strain, metadata$strain)])
  
  counter <- 1e6
  
  for (j in 2:length(unique(df.fragments.cumulative$seqnames))) {
    
    df.fragments.cumulative[j, c("start2", "end2")] <- df.fragments.cumulative[j, c("start2", "end2")] + counter
    
    counter <- counter + 1e6
    
  }
  
  #Get starships
  df.starships <- starship.elements %>%
    filter(strain == strains$strain[i]) %>%
    left_join(df.fragments.cumulative, by="seqnames") %>%
    mutate(plot.xmin=ifelse(start2==1, start.x, start.x+start2),
           plot.xmax=ifelse(end2==1, end.x, end.x+start2))
  
  #Get CSEPs
  CSEPs <- sub("\\..*", "", unlist(strsplit(unlist(as.vector(CSEP.genes[strains$file2[i]])), split=",")))
  
  #Get genes
  df.genes <- as.data.frame(gff) %>%
    filter(type == "gene", biotype != "transposable_element_gene") %>%
    mutate(seqnames=sub(".*v1.1_", "", seqnames)) %>%
    filter(seqnames %in% synteny$contig[synteny$strain == strains$file2[i]]) %>%
    left_join(df.fragments.cumulative, by="seqnames") %>%
    mutate(plot.xmin=ifelse(start2==1, start.x, start.x+start2),
           plot.xmax=ifelse(end2==1, end.x, end.x+start2),
           CSEP=ifelse(sub("-", "", ID) %in% CSEPs, "Y", NA))
  
  #Get TE genes
  df.genes.tes <- as.data.frame(gff) %>%
    filter(type == "gene", biotype == "transposable_element_gene") %>%
    mutate(seqnames=sub(".*v1.1_", "", seqnames)) %>%
    filter(seqnames %in% synteny$contig[synteny$strain == strains$file2[i]]) %>%
    left_join(df.fragments.cumulative, by="seqnames") %>%
    mutate(plot.xmin=ifelse(start2==1, start.x, start.x+start2),
           plot.xmax=ifelse(start2==1, end.x, end.x+start2))
  
  #Summarise gene density for 100,000 bp windows of each fragment
  df.genes.split <- df.genes %>%
    split(.$seqnames)
  
  for (j in 1:length(df.genes.split)) {
    
    tmp <- df.genes.split[[j]] %>%
      mutate(
        window=cut(plot.xmin, dig.lab=8,
                   breaks=seq(
                     df.fragments.cumulative$start2[
                       df.fragments.cumulative$seqnames == 
                         unique(df.genes.split[[j]]$seqnames)
                     ],
                     df.fragments.cumulative$end2[
                       df.fragments.cumulative$seqnames == 
                         unique(df.genes.split[[j]]$seqnames)
                     ], 100000))) %>%
      group_by(window) %>%
      summarise(num=n()) %>%
      mutate(tmp=str_sub(window, 2, -2)) %>% 
      separate(tmp, c("min", "max"), sep=",") %>% 
      mutate_at(c("min", "max"), as.double) %>%
      select(num, max) %>%
      filter(!is.na(max)) %>%
      arrange(max) %>%
      complete(max=full_seq(max, period=100000), fill=list(num=0)) %>%
      mutate(seqnames=unique(df.genes.split[[j]]$seqnames))
    
    assign(paste0("tmp.gene.density.", j), tmp)
    
  }
  
  #Combine results for each fragment
  df.gene.density <- do.call("rbind", mget(ls(pattern="tmp.gene.density."))) %>%
    mutate(strain=strains$strain[i],
           new.strain=metadata$new.strain[match(strain, metadata$strain)])
  
  rm(list=ls(pattern="tmp.gene.density."))
  
  #Get CSEP gene locations
  df.CSEPs <- df.genes %>%
    filter(sub("-", "", ID) %in% CSEPs)
  
  #Summarise CSEP density for 100,000 bp windows of each fragment
  df.CSEPs.split <- df.CSEPs %>%
    split(.$seqnames)
  
  for (j in 1:length(df.CSEPs.split)) {
    
    tmp <- df.CSEPs.split[[j]] %>%
      mutate(
        window=cut(plot.xmin, dig.lab=8,
                   breaks=seq(
                     df.fragments.cumulative$start2[
                       df.fragments.cumulative$seqnames == 
                         unique(df.CSEPs.split[[j]]$seqnames)
                     ],
                     df.fragments.cumulative$end2[
                       df.fragments.cumulative$seqnames == 
                         unique(df.CSEPs.split[[j]]$seqnames)
                     ], 100000))) %>%
      group_by(window) %>%
      summarise(num=n()) %>%
      mutate(tmp=str_sub(window, 2, -2)) %>% 
      separate(tmp, c("min", "max"), sep=",") %>% 
      mutate_at(c("min", "max"), as.double) %>%
      select(num, max) %>%
      filter(!is.na(max)) %>%
      arrange(max) %>%
      complete(max=full_seq(max, period=100000), fill=list(num=0)) %>%
      mutate(seqnames=unique(df.CSEPs.split[[j]]$seqnames))
    
    assign(paste0("tmp.CSEP.density.", j), tmp)
    
  }
  
  #Combine results for each fragment
  df.CSEP.density <- do.call("rbind", mget(ls(pattern="tmp.CSEP.density."))) %>%
    mutate(strain=strains$strain[i],
           new.strain=metadata$new.strain[match(strain, metadata$strain)])
  
  rm(list=ls(pattern="tmp.CSEP.density."))
  
  #Get HCN gene positions
  df.hcn <- df.genes %>%
    filter(sub("-", "", ID) %in% sub("\\..*", "", hcn.genes$gene)) %>%
    mutate(hog=hcn.genes$hog[match(sub("-", "", ID), sub("\\..*", "", hcn.genes$gene))])
  
  #Import TE annotation
  te.gff <- import(
    paste0("S:/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Analysis/",
           strains$file2[i], "/Analysis/eirepeat-1.1.0/output/all_interspersed_repeats.gff3")
  )
  
  #Summarise TE density for 100,000 bp windows of each fragment
  df.te.split <- as.data.frame(te.gff) %>%
    filter(type == "match") %>%
    select(seqnames, start, end, Name) %>%
    mutate(seqnames=sub(".*v1.1_", "", seqnames)) %>%
    filter(seqnames %in% synteny$contig[synteny$strain == strains$file2[i]]) %>%
    left_join(df.fragments.cumulative, by="seqnames") %>%
    mutate(plot.xmin=ifelse(start2==1, start.x, start.x+start2),
           plot.xmax=ifelse(end2==1, end.x, end.x+start2)) %>%
    split(.$seqnames)
  
  for (j in 1:length(df.te.split)) {
    
    tmp <- df.te.split[[j]] %>%
      mutate(
        window=cut(plot.xmin, dig.lab=8,
                   breaks=seq(
                     df.fragments.cumulative$start2[
                       df.fragments.cumulative$seqnames == 
                         unique(df.te.split[[j]]$seqnames)
                     ],
                     df.fragments.cumulative$end2[
                       df.fragments.cumulative$seqnames == 
                         unique(df.te.split[[j]]$seqnames)
                     ], 100000))) %>%
      group_by(window) %>%
      summarise(num=n()) %>%
      mutate(tmp=str_sub(window, 2, -2)) %>% 
      separate(tmp, c("min", "max"), sep=",") %>% 
      mutate_at(c("min", "max"), as.double) %>%
      select(num, max) %>%
      filter(!is.na(max)) %>%
      arrange(max) %>%
      complete(max=full_seq(max, period=100000), fill=list(num=0)) %>%
      mutate(seqnames=unique(df.te.split[[j]]$seqnames))
    
    assign(paste0("tmp.te.density.", j), tmp)
    
  }
  
  #Combine results for each fragment
  df.te.density <- do.call("rbind", mget(ls(pattern="tmp.te.density."))) %>%
    mutate(strain=strains$strain[i],
           new.strain=metadata$new.strain[match(strain, metadata$strain)])
  
  rm(list=ls(pattern="tmp.te.density."))
  
  #Get GC content in 100,000 bp windows
  df.gc <- read.csv(paste0("S://gc/", strains$file2[i], "_gc_100000.tsv"), sep="\t") %>%
    dplyr::rename(seqnames="X.1_usercol", gc="X5_pct_gc", y="X3_usercol") %>%
    mutate(seqnames=sub(".*EI_v1.1_", "", seqnames)) %>%
    filter(seqnames %in% synteny$contig[synteny$strain == strains$file2[i]]) %>%
    select(seqnames, gc, y) %>%
    left_join(df.fragments.cumulative, by="seqnames") %>%
    group_by(seqnames) %>%
    mutate(plot.xmax=ifelse(start2==1, y, y+start2),
           diff=plot.xmax-lag(plot.xmax),
           plot.xmin=ifelse(diff < 100000, lag(plot.xmax), plot.xmax-100000),
           plot.xmin=ifelse(is.na(plot.xmin), start2, plot.xmin),
           diff=ifelse(is.na(diff), 100000, diff)) %>%
    filter(diff == 100000)
  
  #Get composite RIP index in 500 bp windows
  df.rip <- read.csv(paste0("S://RIP/", strains$file1[i], "/", strains$file1[i], "_RIP.bed"),
                     sep="\t", header=FALSE) %>%
    dplyr::rename(seqnames="V1", start="V2", end="V3", CRI="V5") %>%
    mutate(seqnames=sub("l", "", seqnames)) %>%
    select(seqnames, start, end, CRI) %>%
    left_join(df.fragments.cumulative, by="seqnames") %>%
    mutate(plot.xmin=ifelse(start2==1, start.x, start.x+start2),
           plot.xmax=ifelse(start2==1, end.x, end.x+start2))
  
  #Get tandem positions
  df.tandems <- get(paste0("df.tandems.", strains$strain[i])) %>%
    mutate(seqnames=contig) %>%
    left_join(df.fragments.cumulative, by="seqnames") %>%
    mutate(plot.xmin=ifelse(start2==1, start.x, start.x+start2),
           plot.xmax=ifelse(start2==1, end.x, end.x+start2))
  
  #Get contigs which are inverted in GENESPACE plot
  inversions <- synteny %>% 
    filter(!is.na(inversion),
           strain == strains$strain[i])
  
  #Invert positions to match
  for (j in inversions$contig) {
    
    df.gc$plot.xmax[df.gc$seqnames == j] <-
      rev(df.gc$plot.xmax[df.gc$seqnames == j])
    
    df.te.density$max[df.te.density$seqnames == j] <-
      rev(df.te.density$max[df.te.density$seqnames == j])
    
    df.gene.density$max[df.gene.density$seqnames == j] <-
      rev(df.gene.density$max[df.gene.density$seqnames == j])
    
    df.CSEP.density$max[df.CSEP.density$seqnames == j] <-
      rev(df.CSEP.density$max[df.CSEP.density$seqnames == j])
    
    if (j %in% df.tandems$seqnames) {
      
      df.tandems$plot.xmin[df.tandems$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.tandems$plot.xmin[df.tandems$seqnames == j]
      
      df.tandems$plot.xmax[df.tandems$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.tandems$plot.xmax[df.tandems$seqnames == j]
      
    }
    
    if (j %in% df.CSEPs$seqnames) {
      
      df.CSEPs$plot.xmin[df.CSEPs$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.CSEPs$plot.xmin[df.CSEPs$seqnames == j]
      
      
      df.CSEPs$plot.xmax[df.CSEPs$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.CSEPs$plot.xmax[df.CSEPs$seqnames == j]
      
    }
    
    if (j %in% df.hcn$seqnames) {
      
      df.hcn$plot.xmin[df.hcn$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.hcn$plot.xmin[df.hcn$seqnames == j]
      
      
      df.hcn$plot.xmax[df.hcn$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.hcn$plot.xmax[df.hcn$seqnames == j]
      
    }
    
    if (j %in% df.genes$seqnames) {
      
      df.genes$plot.xmin[df.genes$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.genes$plot.xmin[df.genes$seqnames == j]
      
      
      df.genes$plot.xmax[df.genes$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.genes$plot.xmax[df.genes$seqnames == j]
      
      df.genes.tes$plot.xmin[df.genes.tes$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.genes.tes$plot.xmin[df.genes.tes$seqnames == j]
      
      
      df.genes.tes$plot.xmax[df.genes.tes$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.genes.tes$plot.xmax[df.genes.tes$seqnames == j]
      
    }
    
    if (j %in% df.starships$seqnames) {
      
      df.starships$plot.xmin[df.starships$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.starships$plot.xmin[df.starships$seqnames == j]
      
      
      df.starships$plot.xmax[df.starships$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.starships$plot.xmax[df.starships$seqnames == j]
      
    }
    
    if (j %in% df.rip$seqnames) {
      
      df.rip$plot.xmin[df.rip$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.rip$plot.xmin[df.rip$seqnames == j]
      
      
      df.rip$plot.xmax[df.rip$seqnames == j] <-
        df.fragments.cumulative$end2[df.fragments.cumulative$seqnames == j] +
        df.fragments.cumulative$start2[df.fragments.cumulative$seqnames == j] -
        df.rip$plot.xmax[df.rip$seqnames == j]
      
    }
    
  }
  
  #Calculate number of genes with isoforms
  df.isoforms <- as.data.frame(gff) %>% 
    filter(type == "mRNA") %>%
    mutate(gene=sub("\\..*", "", ID)) %>%
    group_by(gene) %>%
    summarise(isoforms=n()) %>%
    mutate(strain=strains$strain[i],
           new.strain=metadata$new.strain[match(strain, metadata$strain)])
  
  #Calculate distances of genes to nearest TE
  te.distances.noncseps <- list()
  te.distances.cseps <- list()
  te.distances.hcn <- list()
  te.distances.nonhcn <- list()
  
  for (j in 1:length(df.genes.split)) {
    
    genes <- df.genes.split[[j]] %>% 
      filter(biotype != "transposable_element_gene") %>%
      mutate(ID=sub("$", ".1", sub("-", "", ID)))
    tes <- df.te.split[[names(df.genes.split)[j]]]
    
    for (k in 1:nrow(genes)) {
      
      start <- genes$start.x[k]
      end <- genes$end.x[k]
      
      if (nrow(tes) == 0) {
        
        distance <- NA
        
      } else if (end < tes$start.x[1]) {
        
        distance <- tes$start.x[1] - end
        
      } else if (start > tes$end.x[nrow(tes)]) {
        
        distance <- start - tes$end.x[nrow(tes)]
        
      } else {
        
        left <- tes$end.x[tes$end.x < start]
        right <- tes$start.x[tes$start.x > end]
        
        if (length(left) == 0 || length(right) == 0) {
          
          distance <- 0
          
        } else {
          
          closest.left <- max(tes$end.x[tes$end.x < start])
          closest.right <- min(tes$start.x[tes$start.x > end])
          
          closer <- which.min(c(start - closest.left,
                                closest.right - end))
          
          if (closer == 1) {
            
            distance <- start - closest.left
            
          } else {
            
            distance <- closest.right - end
            
          }
          
        }
        
      }
      
      if (is.na(genes$CSEP[k])) {
        
        te.distances.noncseps[[names(df.genes.split)[j]]][
          length(te.distances.noncseps[[names(df.genes.split)[j]]])+1] <- distance
        
      } else {
        
        te.distances.cseps[[names(df.genes.split)[j]]][
          length(te.distances.cseps[[names(df.genes.split)[j]]])+1
        ] <- distance
        
      }
      
      if (genes$ID[k] %in% hcn.genes$gene) {
        
        te.distances.hcn[[names(df.genes.split)[j]]][
          length(te.distances.hcn[[names(df.genes.split)[j]]])+1
        ] <- distance
        
      } else {
        
        te.distances.nonhcn[[names(df.genes.split)[j]]][
          length(te.distances.nonhcn[[names(df.genes.split)[j]]])+1] <- distance
        
      }
      
    }
    
  }
  
  print("Proportion of non-CSEPs within 5Kb of TE")
  print(length(which(unlist(te.distances.noncseps) < 5000)) / length(unlist(te.distances.noncseps)))
  print("Proportion of CSEPs within 5Kb of TE")
  print(length(which(unlist(te.distances.cseps) < 5000)) / length(unlist(te.distances.cseps)))
  print("Proportion of non-HCN within 5Kb of TE")
  print(length(which(unlist(te.distances.nonhcn) < 5000)) / length(unlist(te.distances.nonhcn)))
  print("Proportion of HCN within 5Kb of TE")
  print(length(which(unlist(te.distances.hcn) < 5000)) / length(unlist(te.distances.hcn)))
  
  #Combine TE distance results
  df.te.distances <- 
    data.frame(group=rep(c("Non-CSEP", "CSEPs", "HCN", "Non-HCN"),
                         c(length(unlist(te.distances.noncseps)),
                           length(unlist(te.distances.cseps)),
                           length(unlist(te.distances.hcn)),
                           length(unlist(te.distances.nonhcn)))),
               distance=c(unlist(te.distances.noncseps),
                          unlist(te.distances.cseps),
                          unlist(te.distances.hcn),
                          unlist(te.distances.nonhcn)),
               new.strain=metadata$new.strain[match(strains$strain[i], metadata$strain)])
  
  #Format dataframes into GRanges objects
  fragments.gr <- toGRanges(
    df.fragments.cumulative %>%
      mutate(
        chr=synteny$pseudochromosome.abb[
          synteny$strain == strains$file2[i]
        ][
          match(df.fragments.cumulative$seqnames,
                synteny$contig[synteny$strain == strains$file2[i]])
        ],
        invert=ifelse(seqnames %in% inversions$contig, "Y", "N"))
  )
  
  #Assign 10,000 bp telomere regions
  telomeres.up <- as.data.frame(rbind(
    as.data.frame(fragments.gr) %>%
      group_by(chr) %>%
      slice_head() %>%
      filter(invert == "N") %>%
      select(seqnames, start, end, chr),
    as.data.frame(fragments.gr) %>%
      group_by(chr) %>%
      slice_tail() %>%
      filter(invert == "Y") %>%
      select(seqnames, start, end, chr)
  ))
  
  telomeres.down <- as.data.frame(rbind(
    as.data.frame(fragments.gr) %>%
      group_by(chr) %>%
      arrange(descending=FALSE) %>%
      slice_tail() %>%
      filter(invert == "N") %>%
      select(seqnames, start, end, chr),
    as.data.frame(fragments.gr) %>%
      group_by(chr) %>%
      slice_head() %>%
      filter(invert == "Y") %>%
      select(seqnames, start, end, chr)
  ))
  
  telomeres.up.gr <- resize(
    toGRanges(telomeres.up),
    10000, fix="start"
  )
  
  telomeres.down.gr <- resize(
    toGRanges(telomeres.down),
    10000, fix="end"
  )
  
  telomeres.updown.gr <- c(telomeres.down.gr, telomeres.up.gr)
  
  #Perform permutation tests on distance of CSEPs to telomeres and TEs
  
  genes.gr <- toGRanges(df.genes)
  cseps.gr <- toGRanges(df.CSEPs)
  
  te.gr <- toGRanges(
    as.data.frame(te.gff) %>%
      filter(type == "match") %>%
      select(seqnames, start, end, Name) %>%
      mutate(seqnames=sub(".*v1.1_", "", seqnames)) %>%
      filter(seqnames %in% synteny$contig[synteny$strain == strains$file2[i]])
  )
  
  set.seed(2)
  
  csep.perm <- permTest(A=cseps.gr, B=telomeres.updown.gr,
                        randomize.function=resampleRegions, universe=genes.gr,
                        mc.set.seed=FALSE,
                        ntimes=1000,
                        evaluate.function=meanDistance,
                        verbose=TRUE)
  
  set.seed(2)
  
  csep.te.perm <- permTest(A=cseps.gr, B=te.gr,
                           randomize.function=resampleRegions, universe=genes.gr,
                           mc.set.seed=FALSE,
                           ntimes=1000,
                           evaluate.function=meanDistance,
                           verbose=TRUE)
  
  assign(paste0("csep.perm.", strains$strain[i]), csep.perm)
  assign(paste0("csep.te.perm.", strains$strain[i]), csep.te.perm)
  
  #Split genes by chromosome
  df.genes.chrsplit <- df.genes %>%
    mutate(chr=synteny$pseudochromosome.abb[synteny$strain == strains$file2[i]][
      match(seqnames, synteny$contig[synteny$strain == strains$file2[i]])
    ]) %>%
    split(.$chr)
  
  #Split CSEPs by chromosome
  df.CSEPs.chrsplit <- df.CSEPs %>%
    mutate(chr=synteny$pseudochromosome.abb[synteny$strain == strains$file2[i]][
      match(seqnames, synteny$contig[synteny$strain == strains$file2[i]])
    ]) %>%
    split(.$chr)
  
  #Split TEs by chromosome
  df.te.chrsplit <- as.data.frame(te.gff) %>%
    filter(type == "match") %>%
    select(seqnames, start, end, Name) %>%
    mutate(seqnames=sub(".*v1.1_", "", seqnames),
           chr=synteny$pseudochromosome.abb[synteny$strain == strains$file2[i]][
             match(seqnames, synteny$contig[synteny$strain == strains$file2[i]])
           ]) %>%
    filter(seqnames %in% synteny$contig[synteny$strain == strains$file2[i]]) %>%
    split(.$chr)
  
  #Perform permutation tests on distance of CSEPs to telomeres and TEs for each chromosome
  for (chr in unique(fragments.gr$chr)) {
    
    genes.gr <- toGRanges(df.genes.chrsplit[[chr]])
    cseps.gr <- toGRanges(df.CSEPs.chrsplit[[chr]])
    
    telomeres.gr <- telomeres.updown.gr[telomeres.updown.gr$chr == chr]
    
    te.gr <- toGRanges(df.te.chrsplit[[chr]])
    
    set.seed(2)
    
    csep.perm <- permTest(A=cseps.gr, B=telomeres.gr,
                          randomize.function=resampleRegions, universe=genes.gr,
                          mc.set.seed=FALSE,
                          ntimes=1000,
                          evaluate.function=meanDistance,
                          verbose=TRUE)
    
    set.seed(2)
    
    csep.te.perm <- permTest(A=cseps.gr, B=te.gr,
                             randomize.function=resampleRegions, universe=genes.gr,
                             mc.set.seed=FALSE,
                             ntimes=1000,
                             evaluate.function=meanDistance,
                             verbose=TRUE)
    
    assign(paste0("csep.perm.", strains$strain[i], ".", chr), csep.perm)
    assign(paste0("csep.te.perm.", strains$strain[i], ".", chr), csep.te.perm)
    
  }
  
  assign(paste0("df.fragments.cumulative.", strains$strain[i]), df.fragments.cumulative)
  assign(paste0("df.gc.cumulative.", strains$strain[i]), df.gc)
  assign(paste0("df.rip.cumulative.", strains$strain[i]), df.rip)
  assign(paste0("df.tandems.cumulative.", strains$strain[i]), df.tandems)
  assign(paste0("df.starships.cumulative.", strains$strain[i]), df.starships)
  assign(paste0("df.genes.cumulative.", strains$strain[i]), df.genes)
  assign(paste0("df.genes.tes.cumulative.", strains$strain[i]), df.genes.tes)
  assign(paste0("df.gene.density.cumulative.", strains$strain[i]), df.gene.density)
  assign(paste0("df.CSEP.density.cumulative.", strains$strain[i]), df.CSEP.density)
  assign(paste0("df.hcn.cumulative.", strains$strain[i]), df.hcn)
  assign(paste0("df.CSEPs.cumulative.", strains$strain[i]), df.CSEPs)
  assign(paste0("df.te.density.cumulative.", strains$strain[i]), df.te.density)
  assign(paste0("df.te.distances.", strains$strain[i]), df.te.distances)
  assign(paste0("df.isoforms.", strains$strain[i]), df.isoforms)
  
}

############################# COMBINE ALL RESULTS ##############################

all.fragments <- do.call("rbind", mget(ls(pattern="df.fragments.cumulative.")))
all.fragments <- plot_orders(all.fragments) %>%
  mutate(colour=factor(colour,
                       levels=c("pseu_chr1", "pseu_chr2", "pseu_chr3", "pseu_chr4",
                                "pseu_chr5", "pseu_chr6", "mixed")))

all.gc <- do.call("rbind", mget(ls(pattern="df.gc.cumulative.")))
all.gc <- plot_orders(all.gc)

all.rip <- do.call("rbind", mget(ls(pattern="df.rip.cumulative.")))

all.genes <- do.call("rbind", mget(ls(pattern="df.genes.cumulative.")))
all.genes <- plot_orders(all.genes)

all.gene.density <- do.call("rbind", mget(ls(pattern="df.gene.density.cumulative.")))
all.gene.density <- plot_orders(all.gene.density)

all.genes.tes <- do.call("rbind", mget(ls(pattern="df.genes.tes.cumulative.")))

all.te.density <- do.call("rbind", mget(ls(pattern="df.te.density.cumulative.")))
all.te.density <- plot_orders(all.te.density)

all.CSEPs <- do.call("rbind", mget(ls(pattern="df.CSEPs.cumulative.")))
all.CSEPs <- plot_orders(all.CSEPs)

all.CSEP.density <- do.call("rbind", mget(ls(pattern="df.CSEP.density.cumulative.")))

all.hcn <- do.call("rbind", mget(ls(pattern="df.hcn.cumulative.")))
all.hcn <- plot_orders(all.hcn)

all.distances <- do.call("rbind", mget(ls(pattern="df.te.distances."))) 
all.distances <- plot_orders_rev(all.distances)

all.isoforms <- do.call("rbind", mget(ls(pattern="df.isoforms.")))

all.starships <- do.call("rbind", mget(ls(pattern="df.starships.cumulative.")))
all.starships <- plot_orders_rev(all.starships)


################# IDEOGRAM - GC content, gene and TE density  ##################

#Make dataframe to highlight putative centromere regions
centromere.df <- data.frame(
  new.strain=c("Gt-LH10", "Gt-LH10", "Gt-LH10", "Gt-LH10",
               "Gt-4e", "Gt-4e", "Gt-4e", "Gt-4e",
               "Gt-23d", "Gt-23d", "Gt-23d",
               "Gt-8d", "Gt-8d", "Gt-8d",
               "Gt-19d1", "Gt-19d1", "Gt-19d1",
               "Ga-3aA1", "Ga-3aA1", "Ga-3aA1", "Ga-3aA1",
               "Ga-CB1", "Ga-CB1", "Ga-CB1", "Ga-CB1",
               "Gh-2C17", "Gh-2C17", "Gh-2C17", "Gh-2C17",
               "Gh-1B17", "Gh-1B17", "Gh-1B17", "Gh-1B17"),
  pos=c(5808638, 28674328, 41378706, 45760099,
        4900001, 15901577, 38616581, 44875057,
        4900001, 40462406, 44692294,
        15505886, 39033017, 43271416,
        15625778, 38388282, 42737699,
        4600001, 14342662, 38093088, 43268553,
        4600001, 15328520, 40339253, 44629335,
        4400001, 10486465, 21059109,  44892291,
        4400001, 10586456, 25889490, 49648259))

#Set plot order
centromere.df <- plot_orders(centromere.df)

#Plot ideograms with GC content, gene and TE density
gg.ideogram <- ggplot(all.fragments, aes(x=end2)) +
  facet_nested(clade.label+new.strain.label~.,
               space="free", switch="y", scales="free_y",
               labeller=label_parsed,
               nest_line=element_line(),
               strip=strip_nested(size="variable")) +
  geom_rect(aes(xmin=start2-200000, xmax=end2+200000,
                ymin=as.numeric(new.strain.label)-0.2,
                ymax=as.numeric(new.strain.label)+0.2,
                fill=as.factor(colour)),
            colour=NA) +
  geom_segment(data=all.gc,
               aes(x=plot.xmax, xend=plot.xmax,
                   y=as.numeric(new.strain.label)+0.05,
                   yend=as.numeric(new.strain.label)+0.15,
                   colour=gc),
               inherit.aes=FALSE) +
  scale_colour_gradient(name="GC content",
                        limits=c(0, 1),
                        low="black", high="white", na.value="transparent") +
  new_scale_colour() +
  geom_segment(data=all.gene.density,
               aes(x=max, xend=max, 
                   y=as.numeric(new.strain.label)-0.05,
                   yend=as.numeric(new.strain.label)+0.05,
                   colour=num),
               inherit.aes=FALSE) +
  scale_colour_gradient(name="Gene density",
                        low="white", high="black", na.value="transparent") +
  new_scale_colour() +
  geom_segment(data=all.te.density,
               aes(x=max, xend=max,
                   y=as.numeric(new.strain.label)-0.15,
                   yend=as.numeric(new.strain.label)-0.05,
                   colour=num),
               inherit.aes=FALSE) +
  scale_colour_gradient(name="TE density",
                        low="white", high="black", na.value="transparent") +
  geom_rect(data=centromere.df,
            aes(ymin=as.numeric(new.strain.label)+0.25,
                ymax=as.numeric(new.strain.label)-0.25,
                xmin=pos-400000, xmax=pos+400000),
            linetype="22",
            colour="black",
            fill=NA,
            linewidth=0.3,
            inherit.aes=FALSE) +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33",
                             "#44BB99", "#77AADD", "dimgrey")) +
  scale_x_continuous(labels=label_number(accuracy=1,
                                         scale=1e-6,
                                         suffix="Mbp"),
                     expand=c(0, 100)) +
  scale_y_continuous(breaks=seq(1, length((levels(all.fragments$new.strain.label)))),
                     labels=levels(all.fragments$new.strain.label)) +
  theme(legend.position=c(0.96, 0.65),
        legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_text(face="bold", size=6),
        legend.key.size=unit(0.2, "cm"),
        legend.margin=margin(0, 0, 0, 0, unit="pt"),
        strip.background=element_blank(),
        strip.text.y.left=element_text(size=6, face="bold", angle=0),
        strip.clip="off",
        panel.spacing=unit(1, "lines"),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank())

#Write to file
pdf(paste0(dir.comp, "ideogram-", Sys.Date(), ".pdf"), width=7, height=5)
gg.ideogram
dev.off()


################## IDEOGRAM - TE density and CSEP locations  ###################

#Plot TE density and CSEP locations
gg.te <- ggplot(all.fragments, aes(x=end2)) +
  facet_nested(clade.label+new.strain.label~.,
               space="free", switch="y", scales="free_y",
               labeller=label_parsed,
               nest_line=element_line(),
               strip=strip_nested(size="variable")) +
  geom_rect(aes(xmin=start2-200000, xmax=end2+200000,
                ymin=as.numeric(new.strain.label)-0.2,
                ymax=as.numeric(new.strain.label)+0.2,
                fill=as.factor(colour)),
            colour=NA) +
  geom_segment(data=all.te.density,
               aes(x=max, xend=max,
                   y=as.numeric(new.strain.label)-0.14,
                   yend=as.numeric(new.strain.label)+0.14,
                   colour=num),
               inherit.aes=FALSE) +
  geom_rect(data=all.CSEPs,
            aes(xmin=plot.xmin, xmax=plot.xmax,
                ymin=as.numeric(new.strain.label)+0.2,
                ymax=as.numeric(new.strain.label)+0.3),
            linewidth=0.05,
            fill="grey",
            colour="black") +
  geom_text(data=all.fragments %>% 
              filter(strain == "Gh-2C17") %>%
              slice_head(n=1),
            aes(x=54000000, y=as.numeric(new.strain.label)+0.2),
            label="CSEPs",
            fontface="bold",
            size=2) +
  geom_segment(data=all.fragments %>% 
                 filter(strain == "Gh-2C17") %>%
                 slice_head(n=1),
               aes(x=51000000, xend=52600000,
                   y=as.numeric(new.strain.label)+0.25,
                   yend=as.numeric(new.strain.label)+0.2),
               linewidth=0.1) +
  scale_colour_gradient(name="TE density",
                        low="white", high="black", na.value="transparent") +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33", 
                             "#44BB99", "#77AADD", "dimgrey")) +
  guides(colour=guide_colourbar(title.position="top",
                                title.theme=element_text(face="bold", size=6)),
         fill=guide_legend(title.position="top",
                           title.theme=element_text(face="bold", size=6))) +
  scale_x_continuous(labels=label_number(accuracy=1,
                                         scale=1e-6,
                                         suffix="Mbp"),
                     expand=c(0, 100)) +
  theme_minimal() +
  theme(legend.position=c(0.95, 0.65),
        legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.key.size=unit(0.2, "cm"),
        legend.margin=margin(0, 0, 0, 0, unit="pt"),
        strip.background=element_blank(),
        strip.text.y.left=element_text(size=6, face="bold", angle=0),
        strip.clip="off",
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        plot.margin=margin(5.5, 5.5, 5.5, 20),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank())


############################ IDEOGRAM - HCN genes  #############################

#Plot all HCN genes
gg.hcn.base <- ggplot(all.fragments, aes(x=end2)) +
  facet_nested(clade.label+new.strain.label~.,
               space="free", switch="y", scales="free_y",
               labeller=label_parsed,
               nest_line=element_line(),
               strip=strip_nested(size="variable")) +
  geom_rect(aes(xmin=start2-200000, xmax=end2+200000,
                ymin=as.numeric(new.strain.label)-0.2,
                ymax=as.numeric(new.strain.label)+0.2,
                fill=as.factor(colour)),
            colour=NA) +
  geom_rect(aes(xmin=start2, xmax=end2,
                ymin=as.numeric(new.strain.label)-0.14,
                ymax=as.numeric(new.strain.label)+0.14),
            fill="white",
            colour=NA) +
  geom_rect(data=all.hcn,
            aes(xmin=plot.xmin, xmax=plot.xmax,
                ymin=as.numeric(new.strain.label)-0.14,
                ymax=as.numeric(new.strain.label)+0.14),
            linewidth=0.05,
            fill="grey",
            colour="black") +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33",
                             "#44BB99", "#77AADD", "dimgrey")) +
  scale_x_continuous(labels=label_number(accuracy=1,
                                         scale=1e-6,
                                         suffix="Mbp"),
                     expand=c(0, 100)) +
  theme(legend.position="none",
        legend.direction="horizontal",
        legend.text=element_text(size=3),
        legend.key.size=unit(0.2, "cm"),
        legend.margin=margin(0, 0, 0, 0, unit="pt"),
        strip.background=element_blank(),
        strip.text.y.left=element_text(size=6, face="bold", angle=0),
        strip.clip="off",
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank())

# #Optionally add labels for all genes (not very informative)
# # gg.hcn.all <- gg.hcn.base +
#   geom_text_repel(data=all.hcn,
#                   aes(label=as.numeric(sub("N0.HOG", "", hog)),
#                       x=(plot.xmin+plot.xmax)/2,
#                       y=as.numeric(new.strain)+0.2),
#                   nudge_y=0.2,
#                   force=0.025,
#                   direction="x",
#                   size=1,
#                   angle=90,
#                   hjust=0,
#                   box.padding=0.01,
#                   segment.size=0.01,
#                   min.segment.length=0,
#                   max.iter=1e4,
#                   max.overlaps=Inf,
#                   max.time=1)

# #Fix order
all.fragments <- all.fragments %>%
  mutate(new.strain=factor(
    new.strain,
    levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
             "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))
  )
all.hcn <- all.hcn %>%
  mutate(new.strain=factor(
    new.strain,
    levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
             "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))
  )

#For each HCN orthogroup...
for (hcn.hog in sort(unique(all.hcn$hog))) {
  
  #Plot ideograms with location of genes
  gg.hcn.hog <- ggplot(all.fragments, aes(x=end2)) +
    geom_rect(aes(xmin=start2, xmax=end2,
                  ymin=as.numeric(new.strain)-0.14,
                  ymax=as.numeric(new.strain)+0.14),
              fill="white",
              colour="dimgrey",
              linewidth=0.1) +
    geom_rect(data=all.hcn %>% filter(hog == hcn.hog),
              aes(xmin=plot.xmin, xmax=plot.xmax,
                  ymin=as.numeric(new.strain)-0.14,
                  ymax=as.numeric(new.strain)+0.14),
              linewidth=0.05,
              fill="red",
              colour="red") +
    labs(title=hcn.hog) +
    scale_x_continuous(labels=label_number(accuracy=1,
                                           scale=1e-6,
                                           suffix="Mbp")) +
    scale_y_continuous(breaks=1:length(levels(all.hcn$new.strain)),
                       labels=all.hcn$clade[match(levels(all.hcn$new.strain),
                                                  all.hcn$new.strain)]) +
    theme(legend.position="none",
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=6),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          plot.title=element_text(size=6, face="bold"),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_blank())
  
  assign(paste0("gg.hcn.hog.", hcn.hog), gg.hcn.hog)
  
}

#Write to file
pdf(paste0(dir.comp, "hcn_facet-", Sys.Date(), ".pdf"), width=7, height=7)
plot_grid(plotlist=mget(ls(pattern="gg.hcn.hog.")))
dev.off()


#Filter for region of clustered HCN orthologues on Gt-8d
df.hcn.zoom <- `df.hcn.cumulative.Gt-8d` %>%
  filter(colour == "pseu_chr3")

#Format RIP for plotting above ideogram
rip.zoom <- all.rip %>%
  filter(new.strain == "Gt-8d") %>%
  mutate(CRI.plot=CRI*0.02)

#Plot HCN region in Gt-8d
gg.hcn.zoom <- ggplot(`df.fragments.cumulative.Gt-8d`,
                      aes(x=end2, y=new.strain)) +
  geom_rect(aes(xmin=start2-150000, xmax=end2+150000,
                fill=as.factor(colour)),
            ymin=0.8,
            ymax=1.2,
            colour=NA) +
  geom_rect(aes(xmin=start2, xmax=end2),
            ymin=0.86,
            ymax=1.14,
            fill="white",
            colour=NA) +
  geom_rect(data=`df.genes.cumulative.Gt-8d`,
            aes(xmin=plot.xmin, xmax=plot.xmax),
            ymin=0.86,
            ymax=1.14,
            linewidth=0.05,
            fill="black",
            colour=NA) +
  geom_rect(data=`df.hcn.cumulative.Gt-8d`,
            aes(xmin=plot.xmin, xmax=plot.xmax),
            ymin=0.86,
            ymax=1.14,
            linewidth=0.05,
            fill="red",
            colour="red") +
  geom_hline(yintercept=1.3, linewidth=0.1, colour="dimgrey") +
  geom_link2(data=rip.zoom,
             aes(x=plot.xmax, y=1.3+CRI.plot, group=seqnames,
                 colour=after_stat(ifelse(y > 1.3, "positive", "negative"))),
             linewidth=0.1) +
  geom_text_repel(data=df.hcn.zoom,
                  aes(label=as.numeric(sub("N0.HOG", "", hog)),
                      x=(plot.xmin+plot.xmax)/2,
                      y=0.86),
                  nudge_y=-0.1,
                  force=1,
                  direction="x",
                  size=2.5,
                  angle=90,
                  hjust=1,
                  box.padding=0.01,
                  segment.size=0.1,
                  min.segment.length=0,
                  max.iter=1e4,
                  max.overlaps=Inf,
                  max.time=1) +
  annotate("text",
           x=(max(df.hcn.zoom$plot.xmax) + min(df.hcn.zoom$plot.xmax)) / 2,
           y=0.6, size=2.5,
           label=paste0(
             comma(max(df.hcn.zoom$start.x) - min(df.hcn.zoom$start.x)),
             " bp")) +
  annotate("segment",
           x=min(df.hcn.zoom$plot.xmax),
           xend=max(df.hcn.zoom$plot.xmax),
           y=0.65, yend=0.65) +
  annotate("text",
           x=26000000, y=0.7,
           label="Gt-8d pseudochromosome 3",
           fontface="bold", size=3.5) +
  annotate("text",
           x=28800000, y=1.35,
           label="RIP",
           fontface="bold", size=2) +
  annotate("text",
           x=28800000, y=1.25,
           label="no RIP",
           fontface="bold", size=2) +
  coord_cartesian(xlim=c(18700000, 29000000)) +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33",
                             "#44BB99", "#77AADD", "dimgrey")) +
  scale_colour_manual(values=c("grey", "black")) +
  scale_x_continuous(labels=label_number(accuracy=1,
                                         scale=1e-6,
                                         suffix="Mbp"),
                     expand=c(0, 100)) +
  theme(legend.position="none",
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(),
        plot.margin=margin(-15, 5.5, -15, 5.5))

#Compare average RIP in region
zoom.range <- range(
  df.hcn.zoom %>%
    filter(sub("-", "", ID) %in% sub("\\.1", "", hcn.genes$gene)) %>%
    pull(plot.xmax)
)

#Whole pseudochromosome
rip.zoom %>%
  summarise(mean(CRI))

#Within region
rip.zoom %>%
  filter(plot.xmin > zoom.range[1],
         plot.xmax < zoom.range[2]) %>%
  summarise(mean(CRI))


########################### IDEOGRAM - avenacinase #############################

#Get avenacinase gene IDs
df.avenacinase <- 
  read.csv("S:/functional_annotation/005_avenacinase/Gt-3aA1_avenacinase.tsv", sep="\t") %>%
  filter(X96.505 > 80)

#Get avenacinase orthogroup
df.avenacinase.ortho <- all.genes %>%
  mutate(hog=all.ortho.genes$hog[
    match(sub("-", "", all.genes$ID), sub("\\.1$", "", all.ortho.genes$gene))
  ]) %>%
  filter(hog == "N0.HOG0008406")

#Plot ideograms with location of avenacinase gene
gg.avenacinase.pos <- ggplot(all.fragments, aes(x=end2)) +
  facet_nested(clade.label+new.strain.label~.,
               space="free", switch="y", scales="free_y",
               labeller=label_parsed,
               nest_line=element_line(),
               strip=strip_nested(size="variable")) +
  geom_rect(aes(xmin=start2-200000, xmax=end2+200000,
                ymin=as.numeric(new.strain.label)-0.2,
                ymax=as.numeric(new.strain.label)+0.2,
                fill=as.factor(colour)),
            colour=NA) +
  geom_rect(aes(xmin=start2, xmax=end2,
                ymin=as.numeric(new.strain.label)-0.14,
                ymax=as.numeric(new.strain.label)+0.14),
            fill="white",
            colour=NA) +
  geom_rect(data=df.avenacinase.ortho,
            aes(xmin=plot.xmin, xmax=plot.xmax,
                ymin=as.numeric(new.strain.label)-0.14,
                ymax=as.numeric(new.strain.label)+0.14),
            linewidth=0.5,
            fill="black",
            colour="black") +
  geom_label(data=df.avenacinase.ortho %>%
               filter(new.strain == "Gt-LH10") %>%
               mutate(pos=(plot.xmin+plot.xmax)/2),
             aes(x=pos, y=as.numeric(new.strain.label)),
             label="avenacinase",
             label.size=NA,
             fontface="bold",
             fill="black",
             colour="white",
             size=2,
             label.padding=unit(2, "pt"),
             vjust=-1) +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33",
                             "#44BB99", "#77AADD", "dimgrey")) +
  scale_x_continuous(labels=label_number(accuracy=1,
                                         scale=1e-6,
                                         suffix="Mbp"),
                     expand=c(0, 100)) +
  coord_cartesian(clip="off") +
  theme(legend.position="none",
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        strip.background=element_blank(),
        strip.text.y.left=element_text(size=6, face="bold", angle=0),
        strip.clip="off",
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(),
        plot.margin=margin(20, 5.5, 5.5, 20))

#Read in avenacinase gene tree
avenacinase.tree <- 
  read.tree(paste0(dir.genephylo, "raxmlng/gaeumannomyces_avenacinase.raxml.support"))

#Truncate excessively long branch
shortened.edge <- avenacinase.tree$edge[which.max(avenacinase.tree$edge.length), 2]
avenacinase.tree$edge.length[which.max(avenacinase.tree$edge.length)] <-
  avenacinase.tree$edge.length[which.max(avenacinase.tree$edge.length)] / 4

#Plot avenacinase gene tree
gg.tree.avenacinase <- ggtree(avenacinase.tree, layout="daylight") +
  xlim(0, 0.075) +
  ylim(0, 0.06) +
  annotate(geom="label",
           x=0.0193,
           y=0.02,
           label="=",
           size=3,
           label.padding=unit(0, "pt"),
           label.size=0) +
  coord_cartesian(clip="off")

#View default tip label positions
gg.tree.avenacinase$data %>%
  filter(isTip == TRUE)

#Create custom label positions to avoid overlap
avenacinase.labels <- 
  data.frame(label=gg.tree.avenacinase$data %>%
               filter(isTip == TRUE) %>%
               pull(label),
             x=c(0.016, 0.008, 0.0071, 0.002, 0.002,
                 0.008, 0.016, 0.0668, 0.057, 0.065),
             y=c(0.049, 0.0495, 0.04, 0.046, 0.043,
                 0.004, 0.001, 0.026, 0.06, 0.057))

#Add tip labels
gg.tree.avenacinase <- gg.tree.avenacinase +
  geom_text(data=avenacinase.labels, 
            aes(x=x, y=y, label=sub("_EIv1.*", "", label)),
            size=2)

#Read in avenacinase alignment
avenacinase.aln <- fa_read(paste0(dir.genephylo, "alignments/avenacinase_genetree_aln.fasta"))

#Split into chunks for plotting
pos <- c(seq(1, unique(width(avenacinase.aln)), 200), unique(width(avenacinase.aln))+1)

for (i in 1:(length(pos)-1)) {
  
  gg.msa <- ggmsa(avenacinase.aln,
                  pos[i], pos[i+1]-1,
                  ref="AAB09777.1_reference",
                  char_width=0.5, seq_name=TRUE,
                  consensus_views=TRUE, disagreement=TRUE, ) +
    geom_seqlogo() +
    scale_y_discrete(expand=c(0, 0)) +
    scale_x_continuous(breaks=seq(pos[i], pos[i+1]-1)) +
    coord_cartesian(xlim=c(pos[i], pos[i]+199), clip="off") +
    theme(axis.text.x=element_text(size=2, angle=90, vjust=0.5, margin=margin(t=-2)),
          axis.text.y=element_text(size=3, margin=margin(r=-20)),
          plot.margin=margin(10, -10, 3, 5.5))
  
  assign(paste0("gg.msa.", formatC(i, width=2, flag=0)), gg.msa)
  
}

gg.avenacinase.aln <- 
  plot_grid(plotlist=mget(ls(pattern="gg.msa.")), align="v", axis="l", ncol=1)

#Write to file
pdf(paste0(dir.comp, "avenacinase-", Sys.Date(), ".pdf"), width=7, height=6.5)
plot_grid(plot_grid(gg.avenacinase.pos, gg.tree.avenacinase
                    , nrow=1, rel_widths=c(1.1, 1), labels="auto"),
          gg.avenacinase.aln, ncol=1, rel_heights=c(1, 2), labels=c("", "c"))
dev.off()


############################ IDEOGRAM - MAT loci ###############################

#Get MAT idiomorph genes
df.mat.loci <- do.call("rbind", lapply(
  Sys.glob("S:/functional_annotation/005_MAT/*_MAT.tsv"),
  function(fn) 
    data.frame(read.csv(fn, sep="\t", header=FALSE))
)) %>%
  dplyr::rename(accession="V1", ID="V2", evalue="V3",
                bitscore="V4", pident="V5", length="V6") %>%
  mutate(locus=recode(accession,
                      `BAC65091.1`="MAT1-1-1",
                      `BAC65094.1`="MAT1-2-1",
                      `ESU14390.1`="MAT1-1-1"),
         ID=sub("\\..*$", "", ID)) %>%
  dplyr::distinct(locus, ID, .keep_all=TRUE) %>%
  left_join(all.genes, by="ID") 

#Plot ideograms of MAT loci
gg.mat <- ggplot(all.fragments, aes(x=end2)) +
  facet_nested(clade.label+new.strain.label~.,
               space="free", switch="y", scales="free_y",
               labeller=label_parsed,
               nest_line=element_line(),
               strip=strip_nested(size="variable")) +
  geom_rect(aes(xmin=start2-200000, xmax=end2+200000,
                ymin=as.numeric(new.strain.label)-0.2,
                ymax=as.numeric(new.strain.label)+0.2,
                fill=as.factor(colour)),
            colour=NA) +
  geom_rect(aes(xmin=start2, xmax=end2,
                ymin=as.numeric(new.strain.label)-0.14,
                ymax=as.numeric(new.strain.label)+0.14),
            fill="white",
            colour=NA) +
  geom_rect(data=df.mat.loci,
            aes(xmin=plot.xmin, xmax=plot.xmax,
                ymin=as.numeric(new.strain.label)-0.14,
                ymax=as.numeric(new.strain.label)+0.14),
            linewidth=0.5,
            fill="black",
            colour="black") +
  geom_label(data=df.mat.loci,
             aes(x=(plot.xmin+plot.xmax)/2, label=locus),
             y=as.numeric(df.mat.loci$new.strain.label)+0.35,
             label.size=NA,
             fontface="bold",
             fill="black",
             colour="white",
             size=1.2,
             label.padding=unit(1, "pt")) +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33",
                             "#44BB99", "#77AADD", "dimgrey")) +
  scale_x_continuous(labels=label_number(accuracy=1,
                                         scale=1e-6,
                                         suffix="Mbp"),
                     expand=c(0, 100)) +
  coord_cartesian(clip="off") +
  theme(legend.position="none",
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        strip.background=element_blank(),
        strip.text.y.left=element_text(size=6, face="bold", angle=0),
        strip.clip="off",
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(),
        plot.margin=margin(20, 5.5, 5.5, 20))

#Write to file
pdf(paste0(dir.comp, "MAT_loci-", Sys.Date(), ".pdf"), width=3.3, height=2.17)
gg.mat
dev.off()


#################### TE/GENE AND CSEP DENSITY CORRELATION  #####################

#Combine TE and CSEP density by window
df.te.csep <- left_join(all.CSEP.density, all.te.density,
                        by=c("max", "strain", "new.strain", "seqnames")) 

#Set plot order
df.te.csep <- plot_orders_rev(df.te.csep)

#Check for normality
for (strain in unique(df.te.csep$new.strain)) {
  
  plot(
    ggarrange(ggqqplot(df.te.csep$num.x[df.te.csep$new.strain == strain]) +
                ggtitle(paste0(strain, " num.x")),
              ggqqplot(df.te.csep$num.y[df.te.csep$new.strain == strain]) +
                ggtitle(paste0(strain, " num.y")))
  )
  
}

#Calculate Kendall's tau for each strain
df.te.csep.cor <- df.te.csep %>%
  group_by(new.strain, clade) %>%
  cor_test(num.x, num.y, method="kendall") %>%
  filter(p < 0.05)

#Set plot order
df.te.csep.cor <- plot_orders_rev(df.te.csep.cor)

#Function to force axis labels to be integers
integer_breaks <- function(n=5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

#Plot scatterplots of TE and CSEP density
gg.te.cseps <- ggplot(df.te.csep, aes(x=num.x, y=num.y)) +
  facet_nested_wrap(~clade.label+new.strain.label, nrow=1,
                    nest_line=element_line(),
                    labeller=label_parsed) +
  geom_point(size=0.3) +
  geom_smooth(method="loess", linewidth=0.3,
              colour="#87AE88", fill="#87AE88") +
  geom_text(data=df.te.csep.cor,
            aes(label=paste('tau', "==", round(cor, 2))),
            parse=TRUE,
            x=5, y=max(na.omit(df.te.csep$num.y)),
            size=1.7, 
            inherit.aes=FALSE) +
  labs(x="CSEP density", y="TE density") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)),
                     breaks=integer_breaks()) +
  coord_cartesian(ylim=c(0, NA)) +
  theme_minimal() +
  theme(strip.text=element_text(size=8),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.spacing=unit(0.8, "lines"),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=5),
        panel.grid.minor=element_blank(),
        plot.margin=margin(5.5, 5.5, 2, 5.5))

#Combine gene and CSEP density by window
df.gene.csep <- left_join(all.CSEP.density, all.gene.density,
                          by=c("max", "strain", "new.strain", "seqnames"))

#Check for normality
for (strain in unique(df.gene.csep$new.strain)) {
  
  plot(
    ggarrange(ggqqplot(df.gene.csep$num.x[df.gene.csep$new.strain == strain]) +
                ggtitle(paste0(strain, " num.x")),
              ggqqplot(df.gene.csep$num.y[df.gene.csep$new.strain == strain]) +
                ggtitle(paste0(strain, " num.y")))
  )
  
}

#Calculate Kendall's tau for each strain
df.gene.csep.cor <- df.gene.csep %>%
  group_by(new.strain) %>%
  cor_test(num.x, num.y, method="kendall") %>%
  filter(p < 0.05)

#Plot scatterplots of gene and CSEP density
gg.gene.csep <- ggplot(df.gene.csep, aes(x=num.x, y=num.y)) +
  facet_wrap(~new.strain, nrow=1) +
  geom_point(size=0.3) +
  geom_smooth(method="loess", linewidth=0.3,
              colour="#87AE88", fill="#87AE88") +
  geom_text(data=df.gene.csep.cor,
            aes(label=paste('tau', "==", round(cor, 2))),
            parse=TRUE,
            x=5, y=8,
            size=1.7) +
  labs(x="CSEP density", y="Gene density") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)),
                     breaks=integer_breaks()) +
  theme_minimal() +
  theme(strip.text=element_blank(),
        panel.spacing=unit(0.8, "lines"),
        axis.title=element_text(size=7),
        axis.text=element_text(size=5),
        panel.grid.minor=element_blank(),
        plot.margin=margin(2, 5.5, 5.5, 5.5))


############################# TE-GENE DISTANCES  ###############################

#Print mean distances (CSEPs versus other genes versus HCN genes) for each strain
all.distances %>%
  group_by(new.strain, group) %>%
  summarise(mean=mean(distance)) %>%
  print(n=36)

#Print mean distances (all genes) for each strain
all.distances %>%
  group_by(new.strain) %>%
  summarise(mean=mean(distance))

#Check for normality
for (strain in unique(all.distances$new.strain)) {
  
  plot(
    ggarrange(
      ggqqplot(all.distances$distance[
        intersect(which(all.distances$new.strain == strain),
                  which(all.distances$group == "CSEPs"))]) +
        ggtitle(paste0(strain, " CSEPs")),
      ggqqplot(all.distances$distance[
        intersect(which(all.distances$new.strain == strain),
                  which(all.distances$group == "Non-CSEP"))]) +
        ggtitle(paste0(strain, " Non-CSEP")),
      ggqqplot(all.distances$distance[
        intersect(which(all.distances$new.strain == strain),
                  which(all.distances$group == "HCN"))]) +
        ggtitle(paste0(strain, " HCN")),
      ggqqplot(all.distances$distance[
        intersect(which(all.distances$new.strain == strain),
                  which(all.distances$group == "Non-HCN"))]) +
        ggtitle(paste0(strain, " Non-HCN"))
    )
  )
  
}

#Do Wilcoxon rank sum test on distance for CSEPs vs other genes for each strain
distances.cseps.test <- all.distances %>%
  filter(group %in% c("Non-CSEP", "CSEPs")) %>%
  group_by(new.strain) %>%
  wilcox_test(distance ~ group) %>%
  mutate(label=ifelse(p < 0.05, "*", ""))

#Set plot order
distances.cseps.test <- plot_orders_rev(distances.cseps.test)

#Do Wilcoxon rank sum test on distance for HCN vs other genes for each strain
distances.hcn.test <- all.distances %>%
  filter(group %in% c("Non-HCN", "HCN"), clade != "GtA") %>%
  group_by(new.strain) %>%
  wilcox_test(distance ~ group) %>%
  mutate(label=ifelse(p < 0.05, "*", ""))

#Set plot order
distances.hcn.test <- plot_orders_rev(distances.hcn.test)

#Check for normality (regardless of CSEP or not)
for (strain in unique(all.distances$new.strain)) {
  
  plot(
    ggqqplot(all.distances %>%
               filter(new.strain == strain, group %in% c("CSEPs", "Non-CSEP")) %>%
               pull(distance)) +
      ggtitle(strain)
  )
  
}

#Do Games-Howell test on distances 
distances.species.test <- all.distances %>%
  filter(group %in% c("CSEPs", "Non-CSEP")) %>%
  mutate(
    new.strain=factor(sub("-", "_", new.strain),
                      levels=c("Gh_1B17", "Gh_2C17", "Ga_CB1", "Ga_3aA1",
                               "Gt_19d1", "Gt_8d", "Gt_23d", "Gt_4e", "Gt_LH10"))
  ) %>%
  games_howell_test(distance ~ new.strain)

#Make letters labels for significance groups
distances.species.test.labels <- data.frame(
  test=
    multcompLetters(
      setNames(distances.species.test$p.adj,
               paste0(distances.species.test$group1,
                      "-", distances.species.test$group2))
    )$Letters,
  strain=
    names(
      multcompLetters(setNames(distances.species.test$p.adj,
                               paste0(distances.species.test$group1,
                                      "-", distances.species.test$group2))
      )$Letters)
) %>%
  mutate(new.strain=sub("_", "-", strain))

#Set plot order
distances.species.test.labels <- plot_orders_rev(distances.species.test.labels)

#Plot box and violin plots for average TE-gene distance
gg.distances.cseps <- ggplot(all.distances %>%
                               filter(group %in% c("Non-CSEP", "CSEPs")),
                             aes(x=new.strain, y=distance, fill=group)) +
  facet_grid(~clade.label, scales="free_x", space="free", switch="x",
             labeller=label_parsed) +
  geom_violin(position=position_dodge(width=0.8),
              linewidth=0.3,
              alpha=0.5,
              show.legend=FALSE) +
  geom_boxplot(width=0.1, linewidth=0.3, outlier.shape=NA,
               position=position_dodge(width=0.8)) +
  geom_signif(data=distances.cseps.test %>% filter(label == "*"),
              aes(xmin=0.8, xmax=1.2, annotations=label, y_position=160000),
              textsize=5, vjust=0.5, size=0.3,
              manual=TRUE,
              inherit.aes=FALSE) +
  geom_text(data=distances.species.test.labels,
            aes(x=new.strain, y=Inf, label=test),
            vjust=1,
            size=2,
            inherit.aes=FALSE) +
  scale_fill_manual(values=c("#C3D6C3", "#4B854C"),
                    labels=c("CSEPs",
                             "Other genes")) +
  scale_colour_manual(values=c("#C3D6C3", "#4B854C"),
                      labels=c("CSEPs",
                               "Other genes")) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)),
                     limits=c(0, NA),
                     labels=comma) +
  labs(y="Distance to closest TE (bp)") +
  coord_cartesian(clip="off") +
  theme_minimal() +
  theme(strip.placement="outside",
        strip.text=element_text(size=7, face="bold"),
        legend.position="top",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=7, margin=margin(0, 5, 0, 0)),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=7),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=7),
        panel.grid.major.x=element_blank(),
        panel.spacing=unit(1, "lines"),
        plot.margin=margin(10, 5.5, 5.5, 5.5))

#Write to file
pdf(paste0(dir.comp, "te_csep_comp-", Sys.Date(), ".pdf"), width=7, height=6.2)
plot_grid(gg.te,
          plot_grid(gg.te.cseps, gg.gene.csep, gg.distances.cseps,
                    align="v", axis="lr", ncol=1,
                    rel_heights=c(1.1, 0.9, 2),
                    labels=c("b", "", "c")),
          rel_heights=c(1, 1.5),
          ncol=1,
          labels=c("a", ""))
dev.off()


#Plot box and violin plots for average TE-HCN distance
gg.distances.hcn <- 
  ggplot(all.distances %>%
           filter(group %in% c("Non-HCN", "HCN"),
                  clade != "GtA") %>%
           droplevels(),
         aes(x=new.strain, y=distance)) +
  facet_grid(~clade.label, scales="free_x", space="free", switch="x",
             labeller=label_parsed) +
  geom_boxplot(aes(colour=group),
               width=0.3, linewidth=0.3,
               position=position_dodge(width=0.8),
               outlier.size=0.3) +
  geom_signif(data=distances.hcn.test %>%
                filter(label == "*") %>%
                mutate(row=row_number()) %>%
                group_by(clade) %>%
                mutate(id=row_number()),
              aes(xmin=id-0.15, xmax=id+0.15, annotations=label, y_position=40000, group=row),
              textsize=5, vjust=0.5, size=0.3,
              tip_length=0.01,
              manual=TRUE,
              inherit.aes=FALSE) +
  scale_colour_manual(values=c("#C3D6C3", "#4B854C"),
                      labels=c("HCN genes",
                               "Other genes")) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)),
                     limits=c(0, NA),
                     labels=comma) +
  labs(y="Distance to closest TE (bp)") +
  coord_cartesian(ylim=c(0, 40500)) +
  theme_minimal() +
  theme(strip.placement="outside",
        strip.text=element_text(size=7, face="bold"),
        legend.position="top",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=7, margin=margin(0, 5, 0, 0)),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=7),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=7),
        axis.line.y=element_line(colour="grey92",
                                 arrow=grid::arrow(length=unit(0.13, "cm"),
                                                   type="closed",
                                                   ends="last")),
        panel.grid.major.x=element_blank(),
        panel.spacing=unit(1, "lines"))

#Write to file
pdf(paste0(dir.comp, "hcn_duplicates-", Sys.Date(), ".pdf"), width=7, height=5)
plot_grid(gg.duplicates, gg.distances.hcn, gg.hcn.zoom,
          nrow=3, rel_heights=c(1, 1.3, 1), labels="auto")
dev.off()


############################ ISOFORM HISTOGRAMS  ###############################

#Calculate total number of genes
all.gene.numbers <- all.isoforms %>%
  group_by(new.strain) %>%
  summarise(total.genes=n_distinct(gene))

#Calculate proportion of total genes with different numbers of isoforms for each strain
all.isoform.numbers <- all.isoforms %>%
  filter(isoforms > 1) %>%
  group_by(new.strain, isoforms) %>%
  summarise(num=n()) %>%
  left_join(all.gene.numbers) %>%
  mutate(prop=num/total.genes)

#Set plot order
all.isoform.numbers <- plot_orders_rev(all.isoform.numbers)

#Calculate overall proportion of total genes with 2+ isoforms for each strain
all.isoform.rates <- all.isoform.numbers %>%
  group_by(new.strain) %>%
  summarise(prop=sum(prop), x=max(isoforms))

#Set plot order
all.isoform.rates <- plot_orders_rev(all.isoform.rates)

#Function to force axis labels to be multiples of 4
two_breaks <- function(...) {
  fxn <- function(x) {
    breaks <- seq(2, max(x), 4)
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

#Plot histogram of number of isoforms
gg.isoforms <- ggplot(all.isoform.numbers, aes(x=isoforms, y=prop*100)) +
  facet_nested(~clade.label+new.strain.label, scales="free_x", space="free",
               nest_line=element_line(), labeller=label_parsed) +
  geom_bar(stat="identity") +
  geom_label(data=all.isoform.rates,
             aes(x=x, y=20, label=paste0(round(prop*100), "%")),
             position=position_nudge(x=-3, y=-1),
             fontface="bold",
             label.padding=unit(0.1, "lines"),
             label.size=0,
             fill="dimgrey",
             colour="white",
             size=2) +
  labs(x="Number of isoforms", y="Proportion of total genes (%)") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
  scale_x_continuous(breaks=two_breaks()) +
  coord_cartesian(ylim=c(0, NA), clip="off") +
  theme_minimal() +
  theme(strip.placement="outside",
        strip.text=element_text(size=7, face="bold"),
        legend.position="top",
        legend.margin=margin(0, 0, 0, 0),
        legend.key.size=unit(8, "pt"),
        legend.title=element_blank(),
        legend.text=element_text(size=7, margin=margin(0, 5, 0, 0)),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=7),
        axis.text.y=element_text(size=5),
        axis.text.x=element_text(size=7),
        panel.spacing=unit(1, "lines"),
        plot.margin=margin(5.5, 5.5, 0, 5.5))

#Write to file
pdf(paste0(dir.comp, "isoforms-", Sys.Date(), ".pdf"), width=7, height=2)
gg.isoforms
dev.off()


######################## CSEP LOCATION PERMUTATIONS  ###########################

#Overall telomere results
unlist(lapply(mget(paste0("csep.perm.", strains$strain)),
              function(x) x$meanDistance$pval))
unlist(lapply(mget(paste0("csep.perm.", strains$strain))
              , function(x) x$meanDistance$zscore))

#Overall TE results
unlist(lapply(mget(paste0("csep.te.perm.", strains$strain)),
              function(x) x$meanDistance$pval))
unlist(lapply(mget(paste0("csep.te.perm.", strains$strain)),
              function(x) x$meanDistance$zscore))

#Summarise per pseudochromosome permutation results
perm.results <- 
  rbind(
    data.frame(
      test="telomeres",
      strain=sub(".pseu.*", "", sub("csep.perm.", "", ls(pattern="csep.perm.*.pseu.*"))),
      chr=sub(".*.pseu", "pseu", sub("csep.perm.", "", ls(pattern="csep.perm.*.pseu.*"))),
      p=unlist(lapply(mget(ls(pattern="csep.perm.*.pseu.*")),
                      function(x) x$meanDistance$pval)),
      z=unlist(lapply(mget(ls(pattern="csep.perm.*.pseu.*")),
                      function(x) x$meanDistance$zscore)),
      direction=unlist(lapply(mget(ls(pattern="csep.perm.*.pseu.*")),
                              function(x) x$meanDistance$alternative))
    ),
    data.frame(
      test="TEs",
      strain=sub(".pseu.*", "", sub("csep.te.perm.", "", ls(pattern="csep.te.perm.*.pseu.*"))),
      chr=sub(".*.pseu", "pseu", sub("csep.te.perm.", "", ls(pattern="csep.te.perm.*.pseu.*"))),
      p=unlist(lapply(mget(ls(pattern="csep.te.perm.*.pseu.*")),
                      function(x) x$meanDistance$pval)),
      z=unlist(lapply(mget(ls(pattern="csep.te.perm.*.pseu.*")),
                      function(x) x$meanDistance$zscore)),
      direction=unlist(lapply(mget(ls(pattern="csep.te.perm.*.pseu.*")),
                              function(x) x$meanDistance$alternative))
    )
  )

#Format telomere distance results for plotting 
perm.telomeres.df <- perm.results %>%
  filter(test == "telomeres") %>%
  mutate(chr=sub("B", "", chr),
         new.strain=factor(metadata$new.strain[match(strain, metadata$strain)], 
                           levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1",
                                    "Gt-19d1", "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10")),
         label=ifelse(p < 0.001, "<0.001", round(p, digits=3)))

#Format TE distance results for plotting
perm.tes.df <- perm.results %>%
  filter(test == "TEs") %>%
  mutate(chr=sub("B", "", chr),
         new.strain=factor(metadata$new.strain[match(strain, metadata$strain)], 
                           levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1",
                                    "Gt-19d1", "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10")),
         label=ifelse(p < 0.001, "<0.001", round(p, digits=3)))

#Plot grid of p vales for telomere distances
gg.perm.telomeres <- ggplot(perm.telomeres.df, aes(y=new.strain, x=chr)) +
  geom_tile(aes(fill=z, alpha=p<0.05), colour="grey97", linewidth=1) +
  geom_point(aes(colour=p<0.05),
             shape=15, alpha=0) +
  geom_text(aes(label=label, colour=p<0.05),
            size=1.8,
            show.legend=FALSE) +
  scale_fill_gradient2(low="#4B854C", mid="white", high="#E69F00",
                       breaks=c(-6, -4, -2, 0, 2, 4),
                       limits=c(-6, 4),
                       guide=guide_colourbar(title="z-score",
                                             title.position="top")) +
  scale_alpha_manual(values=c(0, 1), guide=NULL) +
  scale_colour_manual(values=c("grey", "black"),
                      guide=guide_legend(title.position="top",
                                         override.aes=list(alpha=1))) +
  scale_x_discrete(position="top",
                   labels=c("chr1/1B", "chr 2/2B",
                            "chr3/3B", "chr4",
                            "chr5/5B", "chr6")) +
  coord_fixed() +
  theme_minimal() +
  theme(legend.margin=margin(0, 10, 0, 0),
        legend.key.size=unit(7, "pt"),
        legend.title=element_text(face="bold", size=6),
        legend.text=element_text(size=5),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(colour="black", size=6),
        axis.text.x=element_text(colour="black", size=5))

#Plot grid of p vales for TE distances
gg.perm.tes <- ggplot(perm.tes.df, aes(y=new.strain, x=chr)) +
  geom_tile(aes(fill=z, alpha=p<0.05), colour="grey97", linewidth=1) +
  geom_point(aes(colour=p<0.05),
             shape=15, alpha=0) +
  geom_text(aes(label=label, colour=p<0.05),
            size=1.8,
            show.legend=FALSE) +
  scale_fill_gradient2(low="#4B854C", mid="white", high="#E69F00",
                       breaks=c(-6, -4, -2, 0, 2, 4),
                       limits=c(-6, 4),
                       guide=guide_colourbar(title="z-score",
                                             title.position="top")) +
  scale_alpha_manual(values=c(0, 1), guide=NULL) +
  scale_colour_manual(values=c("grey", "black"),
                      guide=guide_legend(title.position="top",
                                         override.aes=list(alpha=1))) +
  scale_x_discrete(position="top",
                   labels=c("chr1/1B", "chr 2/2B",
                            "chr3/3B", "chr4",
                            "chr5/5B", "chr6")) +
  coord_fixed() +
  theme_minimal() +
  theme(legend.margin=margin(0, 10, 0, 0),
        legend.key.size=unit(7, "pt"),
        legend.title=element_text(face="bold", size=6),
        legend.text=element_text(size=5),
        panel.grid=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_text(colour="black", size=6),
        axis.text.x=element_text(colour="black", size=5))

#Write to file
pdf(paste0(dir.comp, "permutations-", Sys.Date(), ".pdf"), width=7, height=3.5)
ggarrange(gg.perm.tes, gg.perm.telomeres,
          labels="auto", common.legend=TRUE, legend="right") 
dev.off()


################################# STARSHIPS ####################################

#Filter out 'empty' Starship
all.starships <- all.starships %>%
  filter(elementID != "Gt-4e_s00063")

#Filter for captain genes
all.captains <- all.starships %>% 
  filter(category == "cap") %>%
  mutate(gene=sub("_", "_EIv1_", gene))

#Read in captain genes identified manually
manual.captains <- 
  read.csv(paste0(dir.comp, "captainSummary_DPS_231106.csv"), header=FALSE) %>%
  dplyr::rename(ID="V1", pos="V2") %>%
  left_join(all.genes)

#Set plot order
manual.captains <- plot_orders_rev(manual.captains)

#Read in all tyrosine recombinases predicted by starfish
all.tyr <- 
  read.csv(paste0(dir.comp, "starfish/geneFinder/gaeumannomyces.bed"),
           sep="\t", header=FALSE) %>%
  dplyr::rename(seqnames="V1", start="V2", end="V3", pos="V2",
                gene="V4", category="V5", orientation="V6", tyrID="V7") %>%
  mutate(gene=sub("_", "_EIv1_", gene))


## Euler ##

#Make list for collecting overlapping groups for euler
starship.methods <- 
  c("starfish tyr"=
      length(Reduce(setdiff, 
                    list(all.tyr$gene,
                         all.captains$gene,
                         manual.captains$ID))),
    "manual\ncaps"=
      length(Reduce(setdiff, 
                    list(manual.captains$ID,
                         all.tyr$gene,
                         all.captains$gene))),
    "starfish\ncaps&starfish tyr"=
      length(setdiff(intersect(all.captains$gene, all.tyr$gene),
                     manual.captains$ID)),
    "starfish tyr&manual\ncaps"=
      length(setdiff(intersect(all.tyr$gene, manual.captains$ID),
                     all.captains$gene)),
    "starfish tyr&starfish\ncaps&manual\ncaps"=
      length(Reduce(intersect, 
                    list(manual.captains$ID,
                         all.tyr$gene,
                         all.captains$gene))))

#Plot euler
set.seed(1)
starship.euler <- euler(starship.methods, shape="ellipse")
eulerr_options(padding=unit(2, "pt"))
starship.euler.plot <- 
  plot(starship.euler,
       labels=list(cex=0.5,
                   col="black",
                   padding=unit(c(1, 1), "pt")),
       fills=list(fill=c("white", "#FFAD48", "red", "khaki1",
                         rep("white", 2), "black")),
       edges=list(col="dimgrey",
                  lwd=0.3),
       quantities=list(labels=c(starship.methods[1], 
                                starship.methods[2],
                                "blank",
                                paste0(starship.methods[4], "*"),
                                starship.methods[3],
                                "blank",
                                starship.methods[5]),
                       cex=0.45,
                       col=c(rep("black", 4), "white")))

#Convert to grob
starship.euler.grob <- as.grob(starship.euler.plot)

## Ideogram ##

#Filter for cargo genes
all.cargo <- all.starships %>% 
  filter(str_count(category, "[.]") == 1)

#Plot ideograms with starship locations
gg.starships.ideogram <- ggplot(all.fragments, aes(x=end2)) +
  geom_rect(aes(xmin=start2-200000, xmax=end2+200000,
                ymin=as.numeric(new.strain)-0.2,
                ymax=as.numeric(new.strain)+0.2,
                fill=as.factor(colour)),
            alpha=0.5,
            colour=NA) +
  geom_rect(aes(xmin=start2, xmax=end2,
                ymin=as.numeric(new.strain)-0.14,
                ymax=as.numeric(new.strain)+0.14),
            fill="white",
            colour=NA,
            linewidth=0.1) +
  geom_rect(data=manual.captains %>%
              filter(!ID %in% all.captains$gene,
                     !ID %in% all.tyr$gene),
            aes(xmin=plot.xmin, xmax=plot.xmax,
                ymin=as.numeric(new.strain)-0.14,
                ymax=as.numeric(new.strain)+0.14),
            linewidth=0.3,
            fill="#FFAD48",
            colour="#FFAD48") +
  geom_rect(data=manual.captains %>%
              filter(!ID %in% all.captains$gene,
                     ID %in% all.tyr$gene),
            aes(xmin=plot.xmin, xmax=plot.xmax,
                ymin=as.numeric(new.strain)-0.14,
                ymax=as.numeric(new.strain)+0.14),
            linewidth=0.3,
            fill="khaki1",
            colour="khaki1") +
  geom_text(data=manual.captains %>%
              filter(!ID %in% all.captains$gene,
                     ID %in% all.tyr$gene),
            aes(x=(plot.xmin+plot.xmax)/2,
                y=as.numeric(new.strain)-0.3,
                label="*"),
            fontface="bold",
            colour="black",
            size=2) +
  geom_rect(data=all.cargo,
            aes(xmin=plot.xmin, xmax=plot.xmax,
                ymin=as.numeric(new.strain)-0.14,
                ymax=as.numeric(new.strain)+0.14),
            linewidth=0.05,
            fill="grey",
            colour="grey") +
  geom_rect(data=all.captains %>% filter(gene %in% manual.captains$ID),
            aes(xmin=plot.xmin, xmax=plot.xmax,
                ymin=as.numeric(new.strain)-0.14,
                ymax=as.numeric(new.strain)+0.14),
            linewidth=0.3,
            fill="black",
            colour="black") +
  geom_rect(data=all.captains %>% filter(!gene %in% manual.captains$ID),
            aes(xmin=plot.xmin, xmax=plot.xmax,
                ymin=as.numeric(new.strain)-0.14,
                ymax=as.numeric(new.strain)+0.14),
            linewidth=0.3,
            fill="red",
            colour="red") +
  geom_label_repel(data=all.captains %>%
                     mutate(label=as.numeric(sub(".*_s", "", elementID))),
                   aes(x=(plot.xmin+plot.xmax)/2,
                       y=as.numeric(new.strain)+0.15,
                       label=label),
                   fill="black",
                   label.size=NA,
                   segment.colour="black",
                   size=2.5,
                   colour="white",
                   fontface="bold",
                   nudge_y=0.36,
                   direction="x",
                   label.padding=unit(2, "pt"),
                   point.padding=NA,
                   box.padding=0,
                   segment.size=0.3,
                   min.segment.length=0) +
  coord_cartesian(clip="off") +
  annotation_custom(starship.euler.grob,
                    xmin=48000000, xmax=58000000,
                    ymin=3.5, ymax=6.5) +
  scale_x_continuous(labels=label_number(accuracy=1,
                                         scale=1e-6,
                                         suffix="Mbp"),
                     expand=c(0, 100)) +
  scale_y_discrete(limits=seq(1, 9.5),
                   labels=levels(all.starships$new.strain)) +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33",
                             "#44BB99", "#77AADD", "dimgrey")) +
  guides(fill=guide_legend(title.position="top",
                           title.theme=element_text(face="bold", size=6))) +
  theme_minimal() +
  theme(legend.position=c(0.97, 0.75),
        legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.key.size=unit(0.2, "cm"),
        legend.margin=margin(0, 0, 0, 0, unit="pt"),
        axis.text.y=element_text(size=8),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        strip.background=element_blank(),
        strip.text.y.left=element_text(size=6, face="bold", angle=0, vjust=-0.1),
        strip.clip="off",
        panel.spacing=unit(0.5, "lines"),
        plot.margin=margin(5.5, 30, 5.5, 30),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank())

#Write to file
pdf(paste0(dir.comp, "starships_locations-", Sys.Date(), ".pdf"), width=7, height=4)
gg.starships.ideogram
dev.off()


## Schematics ##

#Combine TE and gene annotations
all.annotations <- bind_rows(all.genes, all.genes.tes)

#Filter for genes
starship.genes <- starship.elements %>%
  filter(category != "insert", category != "flank", elementID != "Gt-4e_s00063") %>%
  mutate(ID=sub("_", "_EIv1_", gene)) %>%
  left_join(all.annotations, by="ID") %>%
  mutate(category=ifelse("transposable_element_gene" == biotype & !is.na(biotype),
                         "transposable_element_gene", category),
         category=sub("\\.", "gene", category),
         clade=factor(metadata$clade[match(sub("_.*", "", elementID), metadata$strain)],
                      levels=c("GtB", "GtA", "Ga", "Gh"))) %>%
  arrange(factor(new.strain, levels(all.fragments$new.strain))) %>%
  select(elementID, feat_id="ID", category, start, end, clade) %>%
  mutate(seq_id=sub("Gt-3aA1_s00046,Gt-3aA1_s00047", "Gt-3aA1_s00046", elementID),
         seq_id=sub("^Gt-3aA1_s00047$", "Gt-3aA1_s00046", seq_id),
         seq_id=str_replace_all(seq_id, c("Gt14LH10"="Gt-LH10",
                                          "Gt-3aA1"="Ga-3aA1",
                                          "Gt-CB1"="Ga-CB1")),
         seq_id=factor(seq_id, levels=rev(unique(seq_id))))

#Filer for flanking repeats
starship.flanks <- starship.elements %>% 
  filter(category == "flank") %>%
  mutate(ID=sub("_", "_EIv1_", gene)) %>%
  left_join(all.annotations, by="ID") %>%
  arrange(factor(new.strain, levels(all.fragments$new.strain))) %>%
  mutate(clade=factor(metadata$clade[match(sub("_.*", "", elementID), metadata$strain)],
                      levels=c("GtB", "GtA", "Ga", "Gh"))) %>%
  select(elementID, feat_id="ID", category, start, end, clade) %>%
  mutate(category=str_match(feat_id, "\\|(.*?)\\|")[,2],
         seq_id=sub("Gt-3aA1_s00046,Gt-3aA1_s00047", "Gt-3aA1_s00046", elementID),
         seq_id=sub("^Gt-3aA1_s00047$", "Gt-3aA1_s00046", seq_id),
         seq_id=str_replace_all(seq_id, c("Gt14LH10"="Gt-LH10",
                                          "Gt-3aA1"="Ga-3aA1",
                                          "Gt-CB1"="Ga-CB1")),
         seq_id=factor(seq_id, levels=rev(unique(seq_id))))

#Filter for nested starship
starship.nested <- starship.elements %>%
  filter(category != "insert", elementID == "Gt-3aA1_s00046,Gt-3aA1_s00047" | elementID == "Gt-3aA1_s00047") %>%
  select(elementID, category, start, end) %>%
  mutate(seq_id=sub("^Gt-3aA1_s00047$", "Gt-3aA1_s00046", elementID),
         seq_id=str_replace_all(seq_id, c("Gt-3aA1"="Ga-3aA1"))) %>%
  filter(seq_id == "Ga-3aA1_s00046") %>%
  group_by(seq_id) %>%
  summarise(start=min(start), end=max(end)) %>%
  mutate(clade=metadata$clade[match(sub("_.*", "", seq_id), metadata$new.strain)],
         seq_id=as.factor(seq_id))

#Combine genes and repeats
starship.feats <- bind_rows(starship.genes, starship.flanks)

#Plot preliminary schematic
gg.starship.schematic.tmp <- gggenomes(genes=starship.feats, feats=starship.nested) %>%
  flip_seqs(where(
    ~.x$seq_id %in%
      c("Gt-LH10_s00089", "Gt-LH10_s00088", "Gt-LH10_s00085", "Gt-LH10_s00079", "Gt-LH10_s00074",
        "Gt-4e_s00053", "Gt-4e_s00058", "Gt-4e_s00062",
        "Gt-23d_s00109", "Gt-23d_s00107", "Gt-23d_s00103", "Gt-23d_s00099",
        "Gt-19d1_s00091",
        "Ga-3aA1_s00044",
        "Ga-CB1_s00036",
        "Gh-1B17_s00001")
  ))

#Pull sequence coordinates
starship.seqs <- gg.starship.schematic.tmp %>% 
  get_seqs() %>%
  mutate(elementID=starship.genes$elementID[match(seq_id, starship.genes$seq_id)],
         seqnames=starship.elements$seqnames[match(elementID, starship.elements$elementID)],
         strain=starship.elements$strain[match(elementID, starship.elements$elementID)],
         new.strain=metadata$new.strain[match(strain, metadata$strain)],
         clade=factor(metadata$clade[match(sub("_.*", "", elementID), metadata$strain)],
                      levels=c("GtB", "GtA", "Ga", "Gh")),
         clade.label=recode_factor(
           clade,
           "GtB"='paste(italic("Gt"), "B")',
           "GtA"='paste(italic("Gt"), "A")',
           "Ga"='paste(italic("Ga"))',
           "Gh"='paste(italic("Gh"))'))

#Add all tracks to schematic
gg.starship.schematic <- gggenomes(seqs=starship.seqs, genes=starship.feats, feats=starship.nested) +
  facet_grid(clade.label~., scales="free",
             space="free", switch="y",
             labeller=label_parsed) +
  geom_feat(data=gg.starship.schematic.tmp %>%
              pull_feats() %>%
              mutate(clade.label=recode_factor(
                clade,
                "Ga"='paste(italic("Ga"))')),
            colour="khaki1", size=5) +
  geom_feat_text(data=gg.starship.schematic.tmp %>%
                   pull_feats() %>%
                   mutate(clade.label=recode_factor(
                     clade,
                     "Ga"='paste(italic("Ga"))')),
                 label="Ga-3aA1_s00047", size=2, nudge_y=-0.8, nudge_x=130000) +
  geom_seq(linewidth=0.3) +
  geom_bin_label(data=gg.starship.schematic.tmp %>%
                   pull_bins() %>%
                   mutate(clade.label=starship.seqs$clade.label[match(bin_id, starship.seqs$bin_id)]),
                 size=2) +
  geom_gene(data=gg.starship.schematic.tmp %>%
              pull_genes() %>%
              mutate(clade.label=recode_factor(
                clade,
                "GtB"='paste(italic("Gt"), "B")',
                "GtA"='paste(italic("Gt"), "A")',
                "Ga"='paste(italic("Ga"))',
                "Gh"='paste(italic("Gh"))')),
            aes(fill=category),
            shape=0,
            colour=NA) +
  scale_fill_manual(values=c("red", "grey", "dimgrey", "#77AADD", "grey", "grey"),
                    breaks=c("cap", "gene", "transposable_element_gene", "tyr"),
                    limits=c("cap", "gene", "transposable_element_gene", "tyr"),
                    labels=c("cap", "gene", "TE", "tyr")) +
  scale_x_continuous(position="top",
                     breaks=pretty_breaks(),
                     labels=label_number(scale=1e-3, suffix=" Kbp")) +
  coord_cartesian(clip="off") +
  theme(legend.position=c(0.85, 0.85),
        legend.text=element_text(size=7),
        legend.key.size=unit(8, "pt"),
        legend.title=element_blank(),
        legend.margin=margin(0,0,-10,0),
        strip.background=element_blank(),
        strip.text.y=element_text(size=8, face="bold", margin=margin(-5, 22, -5, 0)),
        panel.spacing=unit(0, "lines"))

#Get flanking repeats for custom shape plotting
starship.repeats <- gg.starship.schematic.tmp %>% 
  get_seqs() %>%
  left_join(starship.flanks, by="seq_id") %>%
  mutate(seq_id=factor(seq_id, levels=levels(starship.genes$seq_id)),
         repeat.x=ifelse(xend-x > 0,
                         end.y-start.x,
                         end.x-start.y),
         repeat.xend=ifelse(xend-x > 0,
                            start.y-start.x,
                            end.x-end.y),
         repeat.pos=(repeat.x+repeat.xend)/2,
         clade=factor(metadata$clade[match(sub("_.*", "", seq_id), metadata$new.strain)],
                      levels=c("GtB", "GtA", "Ga", "Gh")),
         clade.label=recode_factor(
           clade,
           "GtB"='paste(italic("Gt"), "B")',
           "GtA"='paste(italic("Gt"), "A")',
           "Ga"='paste(italic("Ga"))',
           "Gh"='paste(italic("Gh"))'))

#Add repeats
gg.starship.schematic <- gg.starship.schematic +
  geom_point(data=starship.repeats,
             aes(x=repeat.pos, y=y, shape="DR/TIR"),
             size=1, stroke=0.2) +
  scale_shape_manual(values=8) +
  guides(fill=guide_legend(order=1),
         shape=guide_legend(order=2,
                            override.aes=list(size=1.6, stroke=0.6)),
         colour=guide_legend(order=3,
                             override.aes=list(linewidth=0.6)))

#Get RIP windows for each element
for (new.strain in unique(starship.seqs$new.strain)) {
  
  elements <- starship.seqs$seq_id[starship.seqs$new.strain == new.strain]
  
  for (element in elements) {
    
    starship.seqs.tmp <- 
      starship.seqs[intersect(which(starship.seqs$new.strain == new.strain),
                              which(starship.seqs$seq_id == element)),]
    
    rip.tmp <- all.rip[intersect(which(all.rip$new.strain == new.strain), which(all.rip$seqnames == starship.seqs.tmp$seqnames)),]
    
    rip.min <- starship.seqs.tmp$start
    rip.max <- starship.seqs.tmp$end
    
    if (starship.seqs.tmp$xend-starship.seqs.tmp$x > 0) {
      
      rip.tmp <- rip.tmp %>%
        filter(start.x > starship.seqs.tmp$start,
               end.x < starship.seqs.tmp$end) %>%
        mutate(rip.x=end.x-rip.min,
               rip.xend=start.x-rip.min,
               seq_id=element)
      
    } else {
      
      rip.tmp <- rip.tmp %>%
        filter(start.x > starship.seqs.tmp$start,
               end.x < starship.seqs.tmp$end) %>%
        mutate(rip.x=rip.max-start.x,
               rip.xend=rip.max-end.x,
               seq_id=element)
      
    }
    
    print(element)
    print(rip.tmp %>% summarise(mean(CRI)))
    
    assign(paste0("rip.starship.", element), rip.tmp)
    
  }
  
}

#Combine element RIP windows
starship.rip <- do.call("rbind", mget(ls(pattern="rip.starship."))) %>%
  mutate(seq_id=factor(seq_id, levels=rev(levels(starship.genes$seq_id)))) %>%
  select(seq_id, rip.x, rip.xend, CRI) %>%
  mutate(baseline=as.numeric(seq_id)+0.4,
         y=baseline+(CRI*0.1),
         clade=factor(metadata$clade[match(sub("_.*", "", seq_id), metadata$new.strain)],
                      levels=c("GtB", "GtA", "Ga", "Gh")),
         clade.label=recode_factor(
           clade,
           "GtB"='paste(italic("Gt"), "B")',
           "GtA"='paste(italic("Gt"), "A")',
           "Ga"='paste(italic("Ga"))',
           "Gh"='paste(italic("Gh"))'))

#Plot RIP track above each element
for (i in 1:length(levels(starship.rip$seq_id))) {
  
  gg.starship.schematic <- gg.starship.schematic +
    geom_segment(data=starship.rip %>% filter(seq_id == levels(seq_id)[i]),
                 aes(x=min(rip.x), xend=max(rip.xend),
                     y=baseline, yend=baseline),
                 linewidth=0.1) +
    geom_link2(data=starship.rip %>% filter(seq_id == levels(seq_id)[i]),
               aes(x=(rip.x+rip.xend)/2, y=y,
                   group=seq_id,
                   colour=stage(baseline, y > colour)),
               linewidth=0.1) +
    scale_colour_manual(values=c("black", "lightgrey"),
                        breaks=c("TRUE", "FALSE"),
                        labels=c("RIP", "no RIP"))
  
}


# Tyr/captain tree #

#Read in tree
tyr.tree <- read.tree(paste0(dir.comp, "tyr_genetree/gaeumannomyces_tyr.raxml.support"))
#Root
tyr.tree <- root(tyr.tree,
                 "OBZ73217.1_Grifola_frondosa",
                 edgelabel=TRUE, resolve.root=TRUE)

#Get metadata
tyr.metadata <- read.csv(paste0(dir.comp, "tyr_genetree/metadata_tyr.csv"))

#Format tip labels
tyr.metadata <- tyr.metadata %>%
  mutate(new.label=paste0('italic("', species, '")'),
         new.label=ifelse(grepl("sp\\.", new.label),
                          sub(' sp\\.")', '"), " sp\\."', new.label),
                          new.label),
         new.label=ifelse(gene != "",
                          paste0(new.label, ', " ', gene, '"'),
                          new.label),
         new.label=ifelse(own == "Y",
                          sub("italic", "bolditalic", paste0('bold(paste(', new.label, '))')),
                          paste0('paste(', new.label, ')')),
         type="published cap",
         type=ifelse(starfish == "captain" & manual == "captain", "cap", type),
         type=ifelse(starfish == "captain" & manual != "captain", "starfish cap", type),
         type=ifelse(starfish != "captain" & manual == "captain", "manual cap", type),
         type=ifelse(starfish == "tyr" & manual == "captain", "starfish tyr, manual cap", type),
         type=ifelse(starfish == "tyr" & manual != "captain", "starfish tyr", type))

#Plot tree
gg.tyr.tree <- ggtree(tyr.tree, linetype=NA) %<+% tyr.metadata +
  geom_tree(aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
            linewidth=0.3,
            show.legend=FALSE) +
  xlim(0, 5) +
  scale_colour_manual(values="grey",
                      na.value="black") +
  geom_tiplab(aes(label=new.label),
              size=1,
              offset=0.05,
              parse=TRUE) +
  geom_tippoint(aes(fill=type),
                size=1,
                stroke=0.3,
                shape=21, colour="dimgrey") +
  scale_fill_manual(breaks=c("published cap", "cap", "starfish cap", "manual cap", "starfish tyr, manual cap", "starfish tyr"),
                    values=c("darkgrey", "black", "red", "#FFAD48", "khaki1", "white")) +
  guides(fill=guide_legend(override.aes=list(size=1.3, stroke=0.5))) +
  coord_cartesian(clip="off") +
  annotation_custom(starship.euler.grob, xmin=0.1, xmax=1.8, ymin=95, ymax=130) +
  theme(legend.position=c(0.25, 0.95),
        legend.title=element_blank(),
        legend.text=element_text(size=6, margin=margin(0, 0, 0, -5)),
        legend.key.size=unit(0.2, "cm"),
        legend.margin=margin(0, 0, 0, 0, unit="pt"))

#Write to file
pdf(paste0(dir.comp, "starships-", Sys.Date(), ".pdf"), width=7, height=6)
plot_grid(gg.tyr.tree, gg.starship.schematic, nrow=1, rel_widths=c(3, 4), labels="auto")
dev.off()

################################################################################

save.image(paste0(dir.comp, "plot_ideograms_", Sys.Date(), ".RData"))
