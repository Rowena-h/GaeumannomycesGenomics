library(rtracklayer)
library(ggplot2)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggnewscale)
library(rstatix)
library(multcompView)
library(cowplot)

#Make list with strains and filenames
strains <- list(strain=c("Gt-19d1", "Gt-8d", "Gt-23d", "Gt14LH10", "Gt-4e",
                         "Gt-CB1", "Gt-3aA1", "Gh-2C17", "Gh-1B17"),
                file1=c("Gt17LH19d1", "Gt17LH48d", "Gt23d", "Gt14LH10", "Gt4e",
                        "GtCB1", "PG3aA1", "NZ1292C17", "Gh1B17"),
                file2=c("Gt-19d1", "Gt-8d", "Gt-23d", "Gt14LH10", "Gt-4e",
                        "Gt-CB1", "Gt-3aA1", "Gh-2C17", "Gh-1B17"))

#Read in synteny information inferred from GENESPACE
synteny <- read.csv("R:/GaeumannomycesGenomics/06_synteny/pseudochromosomes.tsv",
                    sep="\t", header=FALSE) %>%
  dplyr::rename(strain="V1", full.contig="V2", pseudochromosome="V3", inversion="V4") %>%
  mutate(contig=sub(".*_", "", full.contig),
         pseudochromosome.abb=sub("-.*", "", pseudochromosome),
         pseudochromosome.abb=ifelse(grepl("B", pseudochromosome.abb), "mixed", pseudochromosome.abb))

#Read in metadata
metadata <- read.csv(paste0("R:/GaeumannomycesGenomics/05_phylogenomics/raxmlng/metadata.csv"))

#Load orthogroup data
load("R:/GaeumannomycesGenomics/07_comparative_genomics/orthogroup-matrices-2023-07-25.RData")

#Remove TEs from orthogroup matrix
orthogroups <- orthogroups %>%
  rownames_to_column(var="orthogroup") %>%
  mutate(biotype=orthogroups.stats$biotype[match(orthogroup, orthogroups.stats$orthogroup)]) %>%
  filter(is.na(biotype)) %>%
  column_to_rownames(var="orthogroup") %>%
  select(-biotype)


## Identify duplicate gene locations ##

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
  mutate(new.strain=metadata$new.strain[match(strain, metadata$strain)],
         new.strain=factor(new.strain, levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1",
                                                "Gt-19d1", "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10")),
         clade=metadata$clade[match(strain, metadata$strain)],
         clade=factor(clade, levels=c("Gh", "Ga", "GtA", "GtB")))
         
levels(df.duplicates.stats$clade) <- c("Gh"=expression(paste(italic("Gh"))),
           "Ga"=expression(paste(italic("Ga"))),
           "GtA"=expression(paste(italic("Gt"), " type A")),
           "GtB"=expression(paste(italic("Gt"), " type B")))

#Plot boxplots of duplicate set size
gg.duplicates <- ggplot(df.duplicates.stats,
                        aes(x=new.strain, y=num, colour=duplicate)) +
  facet_grid(~clade, scales="free_x", space="free", switch="x",
             labeller=label_parsed) +
  geom_boxplot(width=0.3,
               linewidth=0.3,
               outlier.size=0.3,
               position=position_dodge(width=0.5)) +
  labs(x=NULL, y="Duplicate gene set size") +
  scale_colour_manual(values=c("#C3D6C3", "#87AE88", "#4B854C"),
                      labels=c("inter-chromosomal",
                               "intra-chromosomal\n(including tandems)",
                               "mixed")) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.2)),
                     limits=c(0, NA)) +
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


## GC and TE/gene density ideograms ##

#Get CSEP genes
CSEP.genes <- orthogroups[match(orthogroups.stats$orthogroup[!is.na(orthogroups.stats$CSEP)], rownames(orthogroups)),]

#Make dataframe with orthogroup copy-number
copynum.df <- data.frame(strain=colnames(orthogroups.count),
                         as.data.frame(t(orthogroups.count))) %>%
  gather("hog", "num", -strain) %>%
  mutate(biotype=orthogroups.stats$biotype[match(hog, orthogroups.stats$orthogroup)]) %>%
  filter(is.na(biotype))

#Filter for high copy-number (HCN) orthogroups
hcn.df <- copynum.df %>%
  filter(hog %in% unique(copynum.df %>%
                           filter(num > 10) %>%
                           pull(hog))) %>%
  filter(num > 0) %>%
  mutate(clade=metadata$clade[match(strain, metadata$strain)])

#Get all genes for orthogroups
all.genes <- orthogroups %>%
  rownames_to_column("hog") %>%
  mutate(biotype=orthogroups.stats$biotype[match(hog, orthogroups.stats$orthogroup)]) %>%
  filter(is.na(biotype)) %>%
  gather(strain, gene, -hog) %>%
  separate_rows(sep=", ", gene)

#Filter for genes in HCN orthogroups
hcn.genes <- all.genes %>%
  filter(hog %in% unique(hcn.df$hog))

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
    mutate(colour=synteny$pseudochromosome.abb[synteny$strain == strains$file2[i]],
           end2=cumsum(end),
           start2=lag(end2, default=0)+1,
           strain=strains$strain[i],
           new.strain=metadata$new.strain[match(strain, metadata$strain)])
  
  counter <- 1e6
  
  for (j in 2:length(unique(df.fragments.cumulative$seqnames))) {
    
    df.fragments.cumulative[j, c("start2", "end2")] <- df.fragments.cumulative[j, c("start2", "end2")] + counter
    
    counter <- counter + 1e6
    
  }

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

  }
  
  #Calculate distances of genes to nearest TE
  te.distances.genes <- list()
  te.distances.cseps <- list()
  
  for (j in 1:length(df.genes.split)) {
    
    genes <- df.genes.split[[j]] %>% filter(biotype != "transposable_element_gene")
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
        
        te.distances.genes[[names(df.genes.split)[j]]][
          length(te.distances.genes[[names(df.genes.split)[j]]])+1] <- distance
        
      } else {
        
        te.distances.cseps[[names(df.genes.split)[j]]][
          length(te.distances.cseps[[names(df.genes.split)[j]]])+1
        ] <- distance
        
      }
      
    }
    
  }
  
  print("Proportion of genes within 5Kb of TE")
  print(length(which(unlist(te.distances.genes) < 5000)) / length(unlist(te.distances.genes)))
  print("Proportion of CSEPs within 5Kb of TE")
  print(length(which(unlist(te.distances.cseps) < 5000)) / length(unlist(te.distances.cseps)))
  
  #Combine TE distance results
  df.te.distances <- data.frame(group=rep(c("Other genes", "CSEPs"),
                                          c(length(unlist(te.distances.genes)),
                                            length(unlist(te.distances.cseps)))),
                                distance=c(unlist(te.distances.genes),
                                           unlist(te.distances.cseps)),
                                new.strain=metadata$new.strain[match(strains$strain[i], metadata$strain)])
  
  assign(paste0("df.fragments.cumulative.", strains$strain[i]), df.fragments.cumulative)
  assign(paste0("df.gc.cumulative.", strains$strain[i]), df.gc)
  assign(paste0("df.tandems.cumulative.", strains$strain[i]), df.tandems)
  assign(paste0("df.gene.density.cumulative.", strains$strain[i]), df.gene.density)
  assign(paste0("df.CSEP.density.cumulative.", strains$strain[i]), df.CSEP.density)
  assign(paste0("df.hcn.cumulative.", strains$strain[i]), df.hcn)
  assign(paste0("df.CSEPs.cumulative.", strains$strain[i]), df.CSEPs)
  assign(paste0("df.te.density.cumulative.", strains$strain[i]), df.te.density)
  assign(paste0("df.te.distances.", strains$strain[i]), df.te.distances)
  
}

## Combine results for all strains ##

all.fragments <- do.call("rbind", mget(ls(pattern="df.fragments.cumulative.")))
all.fragments$new.strain <- 
  factor(all.fragments$new.strain,
         levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                  "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))
all.fragments$colour <- 
  factor(all.fragments$colour,
         levels=c("pseu_chr1", "pseu_chr2", "pseu_chr3", "pseu_chr4",
                  "pseu_chr5", "pseu_chr6", "mixed"))

all.gc <- do.call("rbind", mget(ls(pattern="df.gc.cumulative.")))
all.gc$new.strain <- 
  factor(all.gc$new.strain,
         levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                  "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))

all.te.density <- do.call("rbind", mget(ls(pattern="df.te.density.cumulative.")))
all.te.density$new.strain <- 
  factor(all.te.density$new.strain,
         levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                  "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))

all.CSEPs <- do.call("rbind", mget(ls(pattern="df.CSEPs.cumulative.")))
all.CSEPs$new.strain <- 
  factor(all.CSEPs$new.strain,
         levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                  "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))

all.gene.density <- do.call("rbind", mget(ls(pattern="df.gene.density.cumulative.")))
all.gene.density$new.strain <- 
  factor(all.gene.density$new.strain,
         levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                  "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))

all.tandems <- do.call("rbind", mget(ls(pattern="df.tandems.cumulative.")))
all.tandems$new.strain <- 
  factor(all.tandems$new.strain,
         levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                  "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))

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

centromere.df$new.strain <- 
  factor(centromere.df$new.strain,
         levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                  "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))

#Plot ideograms with GC content, gene and TE density
gg.ideogram <- ggplot(all.fragments, aes(x=end2, y=new.strain)) +
  geom_rect(aes(xmin=start2-200000, xmax=end2+200000, fill=as.factor(colour)),
            ymin=as.numeric(all.fragments$new.strain)-0.2,
            ymax=as.numeric(all.fragments$new.strain)+0.2,
            colour=NA) +
  geom_segment(data=all.gc,
               aes(x=plot.xmax, xend=plot.xmax, colour=gc),
               y=as.numeric(all.gc$new.strain)+0.05,
               yend=as.numeric(all.gc$new.strain)+0.15,
               inherit.aes=FALSE) +
  scale_colour_gradient(name="GC content",
                        limits=c(0, 1),
                        low="black", high="white", na.value="transparent") +
  new_scale_colour() +
  geom_segment(data=all.gene.density,
               aes(x=max, xend=max, colour=num),
               y=as.numeric(all.gene.density$new.strain)-0.05,
               yend=as.numeric(all.gene.density$new.strain)+0.05,
               inherit.aes=FALSE) +
  scale_colour_gradient(name="Gene density",
                        low="white", high="black", na.value="transparent") +
  new_scale_colour() +
  geom_segment(data=all.te.density,
               aes(x=max, xend=max, colour=num),
               y=as.numeric(all.te.density$new.strain)-0.15,
               yend=as.numeric(all.te.density$new.strain)-0.05,
               inherit.aes=FALSE) +
  scale_colour_gradient(name="TE density",
                        low="white", high="black", na.value="transparent") +
  geom_rect(data=centromere.df,
            aes(ymin=as.numeric(new.strain)+0.25, ymax=as.numeric(new.strain)-0.25,
                xmin=pos-400000, xmax=pos+400000),
            linetype="22",
            colour="black",
            fill=NA,
            linewidth=0.3,
            inherit.aes=FALSE) +
  annotate(geom="segment",
           x=48500000, xend=50500000,
           y=8.1, yend=8.8,
           linewidth=0.2) +
  annotate(geom="segment",
           x=48500000, xend=50500000,
           y=8, yend=7.6,
           linewidth=0.2) +
  annotate(geom="segment",
           x=48500000, xend=50500000,
           y=7.9, yend=6.4,
           linewidth=0.2) +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33", "#44BB99", "#77AADD", "dimgrey")) +
  scale_x_continuous(labels=label_number(accuracy=1,
                                         scale=1e-6,
                                         suffix="Mbp"),
                     expand=c(0, 100)) +
  theme(legend.position=c(0.96, 0.65),
        legend.direction="vertical",
        legend.text=element_text(size=5),
        legend.title=element_text(face="bold", size=6),
        legend.key.size=unit(0.2, "cm"),
        legend.margin=margin(0, 0, 0, 0, unit="pt"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=8),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        #plot.margin=margin(t=20),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank())

#Write to file
pdf(paste0("R://GaeumannomycesGenomics/07_comparative_genomics/ideogram-", Sys.Date(), ".pdf"), width=7, height=5)
gg.ideogram
dev.off()

#Plot TE density and CSEP locations
gg.te <- ggplot(all.fragments, aes(x=end2, y=new.strain)) +
  geom_rect(aes(xmin=start2-200000, xmax=end2+200000, fill=as.factor(colour)),
            ymin=as.numeric(all.fragments$new.strain)-0.2,
            ymax=as.numeric(all.fragments$new.strain)+0.2,
            colour=NA) +
  geom_segment(data=all.te.density,
               aes(x=max, xend=max, colour=num),
               y=as.numeric(all.te.density$new.strain)-0.14,
               yend=as.numeric(all.te.density$new.strain)+0.14,
               inherit.aes=FALSE) +
  geom_rect(data=all.CSEPs,
            aes(xmin=plot.xmin, xmax=plot.xmax),
            linewidth=0.05,
            fill="grey",
            colour="black",
            ymin=as.numeric(all.CSEPs$new.strain)+0.2,
            ymax=as.numeric(all.CSEPs$new.strain)+0.3) +
  annotate(geom="text",
           label="CSEPs",
           x=54000000, y=2,
           size=2,
           fontface="bold") +
  annotate(geom="segment",
           x=52000000, xend=52800000,
           y=1.34, yend=1.75,
           linewidth=0.1) +
  scale_colour_gradient(name="TE density",
                        #limits=c(0, 1),
                        low="white", high="black", na.value="transparent") +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33", "#44BB99", "#77AADD", "dimgrey")) +
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
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=8),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        plot.margin=margin(5.5, 5.5, 5.5, 20),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank())

#Combine all HCN gene positions
all.hcn <- do.call("rbind", mget(ls(pattern="df.hcn.cumulative.")))
all.hcn$new.strain <- 
  factor(all.hcn$new.strain,
         levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                  "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))

#Plot all HCN genes
gg.hcn.base <- ggplot(all.fragments, aes(x=end2, y=new.strain)) +
  geom_rect(aes(xmin=start2-200000, xmax=end2+200000, fill=as.factor(colour)),
            ymin=as.numeric(all.fragments$new.strain)-0.2,
            ymax=as.numeric(all.fragments$new.strain)+0.2,
            colour=NA) +
  geom_rect(aes(xmin=start2, xmax=end2),
            ymin=as.numeric(all.fragments$new.strain)-0.14,
            ymax=as.numeric(all.fragments$new.strain)+0.14,
            fill="white",
            colour=NA) +
  geom_rect(data=all.hcn,
            aes(xmin=plot.xmin, xmax=plot.xmax),
                ymin=as.numeric(all.hcn$new.strain)-0.14,
                ymax=as.numeric(all.hcn$new.strain)+0.14,
            linewidth=0.05,
            fill="grey",
            colour="black") +
  scale_fill_manual(name="Syntenic block",
                    values=c("#FFAABB", "#EE8866", "#EEDD88", "#BBCC33", "#44BB99", "#77AADD", "dimgrey")) +
  scale_x_continuous(labels=label_number(accuracy=1,
                                         scale=1e-6,
                                         suffix="Mbp"),
                     expand=c(0, 100)) +
  theme(legend.position="none",
        legend.direction="horizontal",
        legend.text=element_text(size=3),
        legend.key.size=unit(0.2, "cm"),
        legend.margin=margin(0, 0, 0, 0, unit="pt"),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=8),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank())

# #Optionally add labels for all genes (not very informative)
# library(ggrepel)
# gg.hcn.all <- gg.hcn.base +
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

#Add clade
all.hcn <- all.hcn %>%
  mutate(clade=metadata$clade[match(new.strain, metadata$new.strain)],
         clade=factor(clade, levels=c("Gh", "Ga", "GtA", "GtB")))

#For each HCN orthogroup...
for (hcn.hog in sort(unique(all.hcn$hog))) {
  
  #Plot ideograms with location of genes
  gg.hcn.hog <- ggplot(all.fragments, aes(x=end2, y=new.strain)) +
    geom_rect(aes(xmin=start2, xmax=end2),
              ymin=as.numeric(all.fragments$new.strain)-0.14,
              ymax=as.numeric(all.fragments$new.strain)+0.14,
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
    scale_y_discrete(labels=all.hcn$clade[match(levels(all.hcn$new.strain), all.hcn$new.strain)]) +
    theme(legend.position="none",
          axis.text=element_blank(),
          axis.text.y=element_text(size=6),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          plot.title=element_text(size=6, face="bold"),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_blank())
  
  assign(paste0("gg.hcn.hog.", hcn.hog), gg.hcn.hog)
  
}

#Write to file (just GtB expansions)
pdf(paste0("R://GaeumannomycesGenomics/07_comparative_genomics/hcn_duplicates-", Sys.Date(), ".pdf"), width=7, height=3.5)
plot_grid(gg.duplicates,
          plot_grid(gg.hcn.hog.N0.HOG0000002 + theme(plot.margin=margin(20, 5.5, 5.5, 5.5)),
                    gg.hcn.hog.N0.HOG0000012 + theme(plot.margin=margin(20, 5.5, 5.5, 5.5)),
                    gg.hcn.hog.N0.HOG0000019 + theme(plot.margin=margin(20, 5.5, 5.5, 5.5)),
                    gg.hcn.hog.N0.HOG0000034 + theme(plot.margin=margin(20, 5.5, 5.5, 5.5)),
                    gg.hcn.hog.N0.HOG0000041 + theme(plot.margin=margin(20, 5.5, 5.5, 5.5)),
                    nrow=1),
          nrow=2, rel_heights=c(1.2, 1), labels="auto")
dev.off()

#Write to file
pdf(paste0("R://GaeumannomycesGenomics/07_comparative_genomics/hcn_facet-", Sys.Date(), ".pdf"), width=7, height=7)
plot_grid(plotlist=mget(ls(pattern="gg.hcn.hog.")))
dev.off()


## TE CSEP correlation ##

#Combine CSEP density results for all strains
all.CSEP.density <- do.call("rbind", mget(ls(pattern="df.CSEP.density.cumulative.")))
all.CSEP.density$new.strain <- 
  factor(all.CSEP.density$new.strain,
         levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1", "Gt-19d1",
                  "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10"))

#Combine TE and CSEP density by window
df.te.csep <- left_join(all.CSEP.density, all.te.density, by=c("max", "strain", "new.strain", "seqnames"))

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
  group_by(new.strain) %>%
  cor_test(num.x, num.y, method="kendall") %>%
  filter(p < 0.05)

#Function to force axis labels to be integers
integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

#Plot scatterplots of TE and CSEP density
gg.te.cseps <- ggplot(df.te.csep, aes(x=num.x, y=num.y)) +
  facet_wrap(~new.strain, nrow=1) +
  geom_point(size=0.3) +
  geom_smooth(method = "loess", linewidth=0.3, colour="#87AE88", fill="#87AE88") +
  geom_text(data=df.te.csep.cor,
            aes(label=paste('tau', "==", round(cor, 2))),
            parse=TRUE,
            x=5, y=max(na.omit(df.te.csep$num.y)),
            size=1.7) +
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
df.gene.csep <- left_join(all.CSEP.density, all.gene.density, by=c("max", "strain", "new.strain", "seqnames"))

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
  cor_test(num.x, num.y, method = "kendall") %>%
  filter(p < 0.05)

#Plot scatterplots of gene and CSEP density
gg.gene.csep <- ggplot(df.gene.csep, aes(x=num.x, y=num.y)) +
  facet_wrap(~new.strain, nrow=1) +
  geom_point(size=0.3) +
  geom_smooth(method="loess", linewidth=0.3, colour="#87AE88", fill="#87AE88") +
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


## TE-gene distances ##

#Combine TE-gene distance results for all strains
all.distances <- do.call("rbind", mget(ls(pattern="df.te.distances."))) %>%
  mutate(new.strain=factor(new.strain,
                           levels=c("Gh-1B17", "Gh-2C17", "Ga-CB1", "Ga-3aA1",
                                    "Gt-19d1", "Gt-8d", "Gt-23d", "Gt-4e", "Gt-LH10")),
         clade=metadata$clade[match(new.strain, metadata$new.strain)],
         clade=factor(clade, levels=c("Gh", "Ga", "GtA", "GtB")))

levels(all.distances$clade) <- c("Gh"=expression(paste(italic("Gh"))),
                                 "Ga"=expression(paste(italic("Ga"))),
                                 "GtA"=expression(paste(italic("Gt"), " type A")),
                                 "GtB"=expression(paste(italic("Gt"), " type B")))

#Print mean distances (CSEPs versus other genes) for each strain
all.distances %>%
  group_by(new.strain, group) %>%
  summarise(mean=mean(distance))

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
                  which(all.distances$group == "Other genes"))]) +
        ggtitle(paste0(strain, " other genes"))
    )
  )
  
}

#Do Wilcoxon rank sum test on distance for CSEPs vs other genes for each strain
distances.test <- all.distances %>%
  group_by(new.strain) %>%
  wilcox_test(distance ~ group) %>%
  mutate(label=ifelse(p < 0.05, "*", ""),
         clade=metadata$clade[match(new.strain, metadata$new.strain)],
         clade=factor(clade, levels=c("Gh", "Ga", "GtA", "GtB")))

levels(distances.test$clade) <- c("Gh"=expression(paste(italic("Gh"))),
                                  "Ga"=expression(paste(italic("Ga"))),
                                  "GtA"=expression(paste(italic("Gt"), " type A")),
                                  "GtB"=expression(paste(italic("Gt"), " type B")))

#Check for normality
for (strain in unique(all.distances$new.strain)) {
  
  plot(
    ggqqplot(all.distances$distance[all.distances$new.strain == strain]) +
      ggtitle(strain)
  )
  
}

#Do Games-Howell test on distances 
distances.species.test <- all.distances %>%
  mutate(new.strain=sub("-", "_", new.strain)) %>%
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
  mutate(new.strain=sub("_", "-", strain),
         clade=metadata$clade[match(new.strain, metadata$new.strain)],
         clade=factor(clade, levels=c("Gh", "Ga", "GtA", "GtB")))

levels(distances.species.test.labels$clade) <- 
  c("Gh"=expression(paste(italic("Gh"))),
    "Ga"=expression(paste(italic("Ga"))),
    "GtA"=expression(paste(italic("Gt"), " type A")),
    "GtB"=expression(paste(italic("Gt"), " type B")))


#Plot box and violin plots for average TE-gene distance
gg.distances <- ggplot(all.distances,
                       aes(x=new.strain, y=distance, fill=group)) +
  facet_grid(~clade, scales="free_x", space="free", switch="x",
             labeller=label_parsed) +
  geom_violin(position=position_dodge(width=0.8),
              linewidth=0.3,
              show.legend=FALSE) +
  geom_boxplot(width=0.1, linewidth=0.3, outlier.size=0.3,
               position=position_dodge(width=0.8)) +
  geom_signif(data=distances.test %>% filter(label == "*"),
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
pdf(paste0("R://GaeumannomycesGenomics/07_comparative_genomics/te_csep_comp-", Sys.Date(), ".pdf"), width=7, height=6)

plot_grid(gg.te,
          plot_grid(gg.te.cseps, gg.gene.csep, gg.distances,
                    align="v", axis="lr", ncol=1,
                    rel_heights=c(1, 1, 2),
                    labels=c("b", "", "c")),
          rel_heights=c(1, 1.5),
          ncol=1,
          labels=c("a", ""))

dev.off()
