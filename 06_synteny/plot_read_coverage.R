library(tidyverse)
library(Gviz)
library(rtracklayer)
library(scales)
library(ggpubr)

options(ucscChromosomeNames=FALSE)

#Read in GENESPACE syntenic blocks
blks <- read.csv("S://genespace_gt_230510/results/syntenicBlock_coordinates.csv")

#Filter for regions of interest
blks.region.1 <- blks %>%
  filter(genome1 == "Gt14LH10" & genome2 == "Gt4e" & chr1 == "Gt14LH10_EI_v1.1_ptg000001") %>%
  filter(chr2 == "Gt4e_EI_v1.1_ptg000005" | chr2 == "Gt4e_EI_v1.1_ptg000003")

blks.region.4 <- blks %>%
  filter(genome1 == "Gt14LH10" & genome2 == "Gt4e" & chr1 == "Gt14LH10_EI_v1.1_ptg000004") %>%
  filter(chr2 == "Gt4e_EI_v1.1_ptg000005" | chr2 == "Gt4e_EI_v1.1_ptg000003")

#Create genome plotting track
gt <- GenomeAxisTrack()

#Add mapped reads track
at.1 <- AlignmentsTrack("S://007_readsmap/Gt14LH10/Gt14LH10_aln.sorted.filtered.bam",
                        isPaired = FALSE, chromosome="ptg000001l", stacking="squish",
                        name="Reads", cex=0.5, cex.title=0.5, cex.axis=0.3, min.height=12)

at.4 <- AlignmentsTrack("S://007_readsmap/Gt14LH10/Gt14LH10_aln.sorted.filtered.bam",
                        isPaired = FALSE, chromosome="ptg000004l", stacking="squish",
                        name="Reads", cex=0.5, cex.title=0.5, cex.axis=0.3, min.height=12)

#Read in the header of the gff3 annotation
header <- readLines(
  paste0("S://CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/Gt14LH10/Gt14LH10_EIv1.release.gff3"),
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
  paste0("S://CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/Gt14LH10/Gt14LH10_EIv1.release.gff3")
)

#Filter for genes
genes <- GRanges(paste0(sub(".*_", "", seqnames(annotation)), "l"), ranges(annotation, use.mcols=TRUE), strand=strand(annotation))

#Add syntenic block for each gene
tmp <- as.data.frame(genes) %>%
  filter(type == "gene" & seqnames == "ptg000001l") %>%
  mutate(chr1=sub("l", "", sub("^", "Gt14LH10_EI_v1.1_", seqnames))) %>%
  inner_join(blks.region.1, by=join_by(chr1, within(x$start, x$end, y$startBp1, y$endBp1)))

genes.1 <- genes[match(tmp$ID, genes$ID),]
genes.1$block <- tmp$chr2

tmp <- as.data.frame(genes) %>%
  filter(type == "gene" & seqnames == "ptg000004l") %>%
  mutate(chr1=sub("l", "", sub("^", "Gt14LH10_EI_v1.1_", seqnames))) %>%
  inner_join(blks.region.4, by=join_by(chr1, within(x$start, x$end, y$startBp1, y$endBp1)))

genes.4 <- genes[match(tmp$ID, genes$ID),]
genes.4$block <- tmp$chr2

#Add track for genes coloured by syntenic block
annoT.1 <- AnnotationTrack(genes.1, name="Genes", stacking="dense",
                           cex.title=0.5, cex=0.5, col="transparent",
                           fill="snow3", lime="#BBCC33", blue="#77AADD",
                           pink="#FFAABB", yellow="#EEDD88", orange="#EE8866")
annoT.4 <- AnnotationTrack(genes.4, name="Genes", stacking="dense",
                           cex.title=0.5, cex=0.5, col="transparent",
                           fill="snow3", lime="#BBCC33", blue="#77AADD",
                           pink="#FFAABB", yellow="#EEDD88", orange="#EE8866")

feature(annoT.1) <- 
  sub("Gt4e_EI_v1.1_ptg000004", "lime",
      sub("Gt4e_EI_v1.1_ptg000007", "blue",
          sub("Gt4e_EI_v1.1_ptg000002", "pink",
              sub("Gt4e_EI_v1.1_ptg000005", "yellow",
                  sub("Gt4e_EI_v1.1_ptg000003", "orange", genes.1$block)))))
feature(annoT.4) <- 
  sub("Gt4e_EI_v1.1_ptg000004", "lime",
      sub("Gt4e_EI_v1.1_ptg000007", "blue",
          sub("Gt4e_EI_v1.1_ptg000002", "pink",
              sub("Gt4e_EI_v1.1_ptg000005", "yellow",
                  sub("Gt4e_EI_v1.1_ptg000003", "orange", genes.4$block)))))

#Write to file
pdf(paste0("R://GaeumannomycesGenomics/06_synteny/Gt-LH10_fusion_coverage_1_", Sys.Date(), ".pdf"),
     width=7, height=4)
plotTracks(list(gt, annoT.1, at.1), from=6400000, to=8000000, sizes=c(1, 1, 8), title.width=0.5)
dev.off()


pdf(paste0("R://GaeumannomycesGenomics/06_synteny/Gt-LH10_fusion_coverage_4_", Sys.Date(), ".pdf"),
     width=7, height=4)
plotTracks(list(gt, annoT.4, at.4), from=500000, to=1500000, sizes=c(1, 1, 8), title.width=0.5)
dev.off()


#Read in tidk results
tidk <- read.csv("S://014_tidk/Gt14LH10/Gt14LH10_telomeric_repeat_windows.tsv", sep="\t") %>%
  filter(id %in% c("ptg000001l", "ptg000004l")) %>%
  mutate(id=sub("ptg000004l", "pseu. chr. 3B", sub("ptg000001l", "pseu. chr. 2B", id))) %>%
  mutate(total=forward_repeat_number + reverse_repeat_number)

#Plot telomeric repeats
gg.tidk <- ggplot(tidk, aes(x=window, y=total)) +
  facet_wrap(~id, nrow=2, strip.position="right") +
  geom_line() +
  scale_x_continuous(labels=label_number(scale=1e-6, suffix="Mbp")) +
  labs(y="Number of telomeric repeats") +
  theme_minimal() +
  theme(strip.text=element_text(size=7),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=7),
        axis.text=element_text(size=5),
        panel.spacing=unit(0.8, "lines"),
        panel.grid.minor=element_blank())

pdf(paste0("R://GaeumannomycesGenomics/06_synteny/Gt-LH10_tidk_", Sys.Date(), ".pdf"),
    width=7, height=2)
ggarrange(gg.tidk, labels="b")
dev.off()
