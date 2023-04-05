library(rtracklayer)
library(ggbio)
library(ggplot2)
library(tidyverse)
library(scales)
library(ggpubr)

#Make list with strains and filenames
strains <- list(strain=c("Gt-19d1", "Gt-8d", "Gt-23d", "Gt-LH10", "Gt-4e",
                         "Gt-CB1", "Gt-3aA1", "Gh-2C17", "Gh-1B17"),
                file1=c("Gt17LH19d1", "Gt17LH48d", "Gt23d", "Gt14LH10", "Gt4e",
                        "GtCB1", "PG3aA1", "NZ1292C17", "Gh1B17"),
                file2=c("Gt-19d1", "Gt-8d", "Gt-23d", "Gt14LH10", "Gt-4e",
                        "Gt-CB1", "Gt-3aA1", "Gh-2C17", "Gh-1B17"))

#For each strain...
for (i in 1:length(strains$strain)) {
  
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
  annotation <- import(
    paste0("S://CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/",
           strains$file2[i], "/", strains$file2[i], "_EIv1.release.gff3")
  )
  
  #Add sequence lengths
  seqlengths(annotation) <- 
    df$X4[match(names(seqlengths(annotation)), df$X2)] -
    df$X3[match(names(seqlengths(annotation)), df$X2)]
  
  #Summarise the fragments with annotated genes and format for plotting
  fragments <- GRanges(seqinfo(annotation))
  names(fragments) <- sub(".*_", "", names(fragments))
  df.fragments <- as.data.frame(fragments)
  df.fragments$seqnames <- sub(".*_", "", df.fragments$seqnames)
  df.fragments$seqnames <- factor(
    df.fragments$seqnames,
    levels=df.fragments$seqnames[order(df.fragments$width, decreasing=FALSE)]
  )
  
  #Read in Tapestry telomere identification results
  telomeres <- read.csv(
    paste0("S://015_tapestry/", strains$file1[i], "/", "contig_details.tsv"),
    sep="\t"
  )
  
  #Format and match to annotation above
  telomeres$Contig <- substr(telomeres$Contig, 1, 9)
  telomeres <- telomeres[match(df.fragments$seqnames, telomeres$Contig),]
  df.telomeres <- telomeres %>%
    gather(telomere, num.telomeres, StartTelomeres, EndTelomeres)
  df.telomeres$num.telomeres[df.telomeres$num.telomeres == 0] <- NA
  
  #Add positions for start and end of fragments
  df.telomeres$x.pos <- match(df.telomeres$Contig, levels(df.fragments$seqnames))
  df.telomeres$y.pos <- ifelse(df.telomeres$telomere == "StartTelomeres", 1, df.telomeres$Length)
  
  #Plot karyograms with telomeres
  gg.karyogram <- ggplot(df.fragments, aes(x=seqnames, y=end)) +
    geom_rect(aes(ymin=start, ymax=end),
              xmin=as.numeric(df.fragments$seqnames)-0.2,
              xmax=as.numeric(df.fragments$seqnames)+0.2,
              fill="white",
              colour="dimgrey") +
    geom_segment(data=df.telomeres,
                 aes(x=x.pos-0.3,
                     xend=x.pos+0.3,
                     y=y.pos,
                     yend=y.pos,
                     colour=num.telomeres),
                 size=0.5) +
    scale_y_continuous(labels=label_number(accuracy=1,
                                           scale=1e-6,
                                           suffix="Mbp"),
                       expand=c(0, 100)) +
    scale_colour_gradient(name="Number of telomeric repeats",
                          limits=c(1, max(na.omit(df.telomeres$num.telomeres))),
                          low="#ffbdbd", high="#ff0000", na.value="transparent") +
    guides(colour=guide_colourbar(title.position="top",
                                  title.theme=element_text(face="bold", size=4))) +
    ggtitle(label=strains$strain[i]) +
    coord_flip(clip="off") +
    theme(legend.position=c(0.8, 1.1),
          legend.direction="horizontal",
          legend.text=element_text(size=3),
          legend.key.size=unit(0.2, "cm"),
          legend.margin=margin(0, 0, 0, 0, unit="pt"),
          axis.text.y=element_text(colour="black",
                                   size=5,
                                   margin=margin(r=5)),
          axis.text.x=element_text(size=5),
          axis.ticks.y=element_blank(),
          axis.title=element_blank(),
          axis.line.x=element_line(),
          plot.title=element_text(size=10, face="bold"),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_blank())
  
  assign(paste0("gg.karyogram.", strains$strain[i]), gg.karyogram)
  
}

#Write to file
tiff(paste0("telomeres_Gt-", Sys.Date(), ".tiff"), width=7, height=4, units="in", compression="lzw", res=600)
ggarrange(plotlist=mget(ls(pattern="gg.karyogram.")))
dev.off()
