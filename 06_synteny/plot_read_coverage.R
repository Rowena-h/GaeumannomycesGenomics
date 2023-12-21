library(tidyverse)
library(Gviz)
library(rtracklayer)
library(scales)
library(ggpubr)

## READ COVERAGE ACROSS TRANSLOCATIONS ##

options(ucscChromosomeNames=FALSE)

#Read in GENESPACE syntenic blocks
blks <- read.csv("S://genespace_gaeumannomyces_230510/results/syntenicBlock_coordinates.csv")

#Filter for regions of interest
blks.Gt14LH10.2B <- blks %>%
  filter(genome1 == "Gt14LH10" & genome2 == "Gt4e" & chr1 == "Gt14LH10_EI_v1.1_ptg000001") %>%
  filter(chr2 == "Gt4e_EI_v1.1_ptg000005" | chr2 == "Gt4e_EI_v1.1_ptg000003")

blks.Gt14LH10.3B <- blks %>%
  filter(genome1 == "Gt14LH10" & genome2 == "Gt4e" & chr1 == "Gt14LH10_EI_v1.1_ptg000004") %>%
  filter(chr2 == "Gt4e_EI_v1.1_ptg000005" | chr2 == "Gt4e_EI_v1.1_ptg000003")

blks.Gh1B17.5B <- blks %>%
  filter(genome1 == "Gh1B17" & genome2 == "GtCB1" & chr1 == "Gh1B17_EI_v1.1_ptg000001") %>%
  filter(chr2 == "GtCB1_EI_v1.1_ptg000001" | chr2 == "GtCB1_EI_v1.1_ptg000003")

blks.Gh1B17.1B <- blks %>%
  filter(genome1 == "Gh1B17" & genome2 == "GtCB1" & chr1 == "Gh1B17_EI_v1.1_ptg000004") %>%
  filter(chr2 == "GtCB1_EI_v1.1_ptg000001" | chr2 == "GtCB1_EI_v1.1_ptg000002" | chr2 == "GtCB1_EI_v1.1_ptg000005")

#Create genome plotting tracks
gt.Gt14LH10 <- GenomeAxisTrack()
gt.Gh1B17 <- GenomeAxisTrack()

#Add mapped reads tracks
at.Gt14LH10.2B <- AlignmentsTrack("S://007_readsmap/Gt14LH10/Gt14LH10_aln.sorted.filtered.bam",
                                  isPaired = FALSE, chromosome="ptg000001l", stacking="squish",
                                  name="Reads", cex=0.5, cex.title=0.5, cex.axis=0.3, min.height=12)

at.Gt14LH10.3B <- AlignmentsTrack("S://007_readsmap/Gt14LH10/Gt14LH10_aln.sorted.filtered.bam",
                                  isPaired = FALSE, chromosome="ptg000004l", stacking="squish",
                                  name="Reads", cex=0.5, cex.title=0.5, cex.axis=0.3, min.height=12)

at.Gh1B17.5B <- AlignmentsTrack("S://007_readsmap/Gh1B17/Gh1B17_aln.sorted.filtered.bam",
                                isPaired = FALSE, chromosome="ptg000001l", stacking="squish",
                                name="Reads", cex=0.5, cex.title=0.5, cex.axis=0.3, min.height=12)

at.Gh1B17.1B <- AlignmentsTrack("S://007_readsmap/Gh1B17/Gh1B17_aln.sorted.filtered.bam",
                                isPaired = FALSE, chromosome="ptg000004l", stacking="squish",
                                name="Reads", cex=0.5, cex.title=0.5, cex.axis=0.3, min.height=12)

#Import the whole annotations
annotation.Gt14LH10 <- import(
  paste0("S://CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/Gt14LH10/Gt14LH10_EIv1.release.gff3")
)

annotation.Gh1B17 <- import(
  paste0("S://CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/Gh-1B17/Gh-1B17_EIv1.release.gff3")
)

#Filter for genes
genes.Gt14LH10 <- GRanges(paste0(sub(".*_", "", seqnames(annotation.Gt14LH10)), "l"),
                          ranges(annotation.Gt14LH10, use.mcols=TRUE),
                          strand=strand(annotation.Gt14LH10))

genes.Gh1B17 <- GRanges(paste0(sub(".*_", "", seqnames(annotation.Gh1B17)), "l"),
                        ranges(annotation.Gh1B17, use.mcols=TRUE),
                        strand=strand(annotation.Gh1B17))

#Add syntenic block for each gene
tmp <- as.data.frame(genes.Gt14LH10) %>%
  filter(type == "gene" & seqnames == "ptg000001l") %>%
  mutate(chr1=sub("l", "", sub("^", "Gt14LH10_EI_v1.1_", seqnames))) %>%
  inner_join(blks.Gt14LH10.2B, by=join_by(chr1, within(x$start, x$end, y$startBp1, y$endBp1)))

genes.Gt14LH10.2B <- genes.Gt14LH10[match(tmp$ID, genes.Gt14LH10$ID),]
genes.Gt14LH10.2B$block <- tmp$chr2

tmp <- as.data.frame(genes.Gt14LH10) %>%
  filter(type == "gene" & seqnames == "ptg000004l") %>%
  mutate(chr1=sub("l", "", sub("^", "Gt14LH10_EI_v1.1_", seqnames))) %>%
  inner_join(blks.Gt14LH10.3B, by=join_by(chr1, within(x$start, x$end, y$startBp1, y$endBp1)))

genes.Gt14LH10.3B <- genes.Gt14LH10[match(tmp$ID, genes.Gt14LH10$ID),]
genes.Gt14LH10.3B$block <- tmp$chr2

tmp <- as.data.frame(genes.Gh1B17) %>%
  filter(type == "gene" & seqnames == "ptg000001l") %>%
  mutate(chr1=sub("l", "", sub("^", "Gh1B17_EI_v1.1_", seqnames))) %>%
  inner_join(blks.Gh1B17.5B, by=join_by(chr1, within(x$start, x$end, y$startBp1, y$endBp1)))

genes.Gh1B17.5B <- genes.Gh1B17[match(tmp$ID, genes.Gh1B17$ID),]
genes.Gh1B17.5B$block <- tmp$chr2

tmp <- as.data.frame(genes.Gh1B17) %>%
  filter(type == "gene" & seqnames == "ptg000004l") %>%
  mutate(chr1=sub("l", "", sub("^", "Gh1B17_EI_v1.1_", seqnames))) %>%
  inner_join(blks.Gh1B17.1B, by=join_by(chr1, within(x$start, x$end, y$startBp1, y$endBp1)))

genes.Gh1B17.1B <- genes.Gh1B17[match(tmp$ID, genes.Gh1B17$ID),]
genes.Gh1B17.1B$block <- tmp$chr2

#Add track for genes coloured by syntenic block
annoT.Gt14LH10.2B <- AnnotationTrack(genes.Gt14LH10.2B, name="Genes", stacking="dense",
                                     cex.title=0.5, cex=0.5, col="transparent",
                                     fill="snow3", yellow="#EEDD88", orange="#EE8866")
annoT.Gt14LH10.3B <- AnnotationTrack(genes.Gt14LH10.3B, name="Genes", stacking="dense",
                                     cex.title=0.5, cex=0.5, col="transparent",
                                     fill="snow3", yellow="#EEDD88", orange="#EE8866")
annoT.Gh1B17.5B <- AnnotationTrack(genes.Gh1B17.5B, name="Genes", stacking="dense",
                                   cex.title=0.5, cex=0.5, col="transparent",
                                   fill="snow3", teal="#44BB99", pink="#FFAABB")
annoT.Gh1B17.1B <- AnnotationTrack(genes.Gh1B17.1B, name="Genes", stacking="dense",
                                   cex.title=0.5, cex=0.5, col="transparent",
                                   fill="snow3", teal="#44BB99", pink="#FFAABB",
                                   orange="#EE8866")

feature(annoT.Gt14LH10.2B) <- 
  sub("Gt4e_EI_v1.1_ptg000005", "yellow", #chr3
      sub("Gt4e_EI_v1.1_ptg000003", "orange", #ch2
          genes.Gt14LH10.2B$block))
feature(annoT.Gt14LH10.3B) <- 
  sub("Gt4e_EI_v1.1_ptg000005", "yellow", #chr3
      sub("Gt4e_EI_v1.1_ptg000003", "orange", #chr2
          genes.Gt14LH10.3B$block))
feature(annoT.Gh1B17.5B) <- 
  sub("GtCB1_EI_v1.1_ptg000001", "teal", #chr5
      sub("GtCB1_EI_v1.1_ptg000003", "pink", #chr1
          genes.Gh1B17.5B$block))
feature(annoT.Gh1B17.1B) <- 
  sub("GtCB1_EI_v1.1_ptg000001", "teal", #chr5
      sub("GtCB1_EI_v1.1_ptg000002", "pink", #chr1
          sub("GtCB1_EI_v1.1_ptg000005", "orange", #chr2
              genes.Gh1B17.1B$block)))

#Write to file
pdf(paste0("R://GaeumannomycesGenomics/06_synteny/Gt-LH10_translocation_coverage_2B_", Sys.Date(), ".pdf"),
    width=7, height=4)
plotTracks(list(gt.Gt14LH10, annoT.Gt14LH10.2B, at.Gt14LH10.2B), from=6400000, to=8000000, sizes=c(1, 1, 8), title.width=0.5)
dev.off()


pdf(paste0("R://GaeumannomycesGenomics/06_synteny/Gt-LH10_translocation_coverage_3B_", Sys.Date(), ".pdf"),
    width=7, height=4)
plotTracks(list(gt.Gt14LH10, annoT.Gt14LH10.3B, at.Gt14LH10.3B), from=500000, to=1500000, sizes=c(1, 1, 8), title.width=0.5)
dev.off()

pdf(paste0("R://GaeumannomycesGenomics/06_synteny/Gh-1B17_translocation_coverage_5B_", Sys.Date(), ".pdf"),
    width=7, height=4)
plotTracks(list(gt.Gh1B17, annoT.Gh1B17.5B, at.Gh1B17.5B), from=3100000, to=4100000, sizes=c(1, 1, 8), title.width=0.5)
dev.off()

pdf(paste0("R://GaeumannomycesGenomics/06_synteny/Gh-1B17_translocation_coverage_1B.1_", Sys.Date(), ".pdf"),
    width=7, height=4)
plotTracks(list(gt.Gh1B17, annoT.Gh1B17.1B, at.Gh1B17.1B), from=400000, to=600000, sizes=c(1, 1, 8), title.width=0.5)
dev.off()

pdf(paste0("R://GaeumannomycesGenomics/06_synteny/Gh-1B17_translocation_coverage_1B.2_", Sys.Date(), ".pdf"),
    width=7, height=4)
plotTracks(list(gt.Gh1B17, annoT.Gh1B17.1B, at.Gh1B17.1B), from=4450000, to=4500000, sizes=c(1, 1, 8), title.width=0.5)
dev.off()


## TIDK TELOMERE IDENTIFICATION ##

#Read in tidk results
tidk.Gt14LH10 <- read.csv("S://014_tidk/Gt14LH10/Gt14LH10_telomeric_repeat_windows.tsv", sep="\t") %>%
  filter(id %in% c("ptg000001l", "ptg000004l")) %>%
  mutate(strain="Gt-LH10",
         id=sub("ptg000004l", "pseu. chr. 3B", sub("ptg000001l", "pseu. chr. 2B", id)),
         total=forward_repeat_number + reverse_repeat_number)

tidk.Gh1B17 <- read.csv("S://014_tidk/Gh1B17/Gh1B17_telomeric_repeat_windows.tsv", sep="\t") %>%
  filter(id %in% c("ptg000001l", "ptg000004l")) %>%
  mutate(strain="Gh-1B17",
         id=sub("ptg000004l", "pseu. chr. 1B", sub("ptg000001l", "pseu. chr. 5B", id)),
         total=forward_repeat_number + reverse_repeat_number)

tidk <- bind_rows(tidk.Gh1B17, tidk.Gt14LH10)

#Plot telomeric repeats
gg.tidk <- ggplot(tidk, aes(x=window, y=total)) +
  facet_wrap(~strain+id, nrow=4, strip.position="right") +
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

pdf(paste0("R://GaeumannomycesGenomics/06_synteny/translocations_tidk_", Sys.Date(), ".pdf"),
    width=7, height=4)
ggarrange(gg.tidk, labels="c")
dev.off()
