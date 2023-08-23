library(GENESPACE)
library(data.table)
library(ggplot2)
library(ape)
library(ggtree)
library(aplot)

## Species tree ##

#Assign outgroup
outgroups <- c("GCA_000193285.1_Mag_poae_ATCC_64411_V1_protein.faa")

#Read in tree metadata
metadata <- read.csv(paste0("../05_phylogenomics/raxmlng/metadata.csv"))
#Format tip labels
metadata$new.label <- ifelse(
  metadata$own == "Y",
  paste0('paste(bolditalic("', metadata$name, '"), bold(" ', metadata$new.strain, '"))'),
  paste0('paste(italic("', metadata$name, '"), " ', metadata$new.strain, '")')
)

#Read in tree
tree <- read.tree("../05_phylogenomics/raxmlng/gaeumannomyces_concat.raxml.support")

#Root tree and remove outgroup
tree <- root(tree, outgroups, edgelabel=TRUE, resolve.root=TRUE)
tree <- drop.tip(tree, outgroups)

#Remove other Gt tip
tree <- drop.tip(tree, "GCF_000145635.1_Gae_graminis_V2_protein.faa_XP")

#Make dataframe of Gt type nodes
types.df <- data.frame(
  type=c("A", "B"),
  node=NA
)

#Find the most recent common ancestor for each type
for (i in 1:length(types.df$type)) {
  
  types.df$node[i] <- MRCA(tree,
                           metadata$tip[metadata$type == types.df$type[i] &
                                          metadata$own == "Y"])
  
}

#Plot  species tree
gg.tree <- ggtree(tree, branch.length="none") %<+% metadata +
  xlim(0, 30) +
  geom_tiplab(aes(label=new.label),
              offset=0.5,
              size=2,
              parse=TRUE)

gg.tree <- gg.tree %>%
  ggtree::rotate(16) %>% ggtree::rotate(17) +
  geom_cladelab(data=types.df,
                mapping=aes(node=node, label=type),
                fontsize=2.5,
                fontface="bold",
                barsize=0.8,
                offset=20,
                offset.text=1)

## Adapted GENESPACE plot function ##

#Make list of strains and filenames
strains <- list(strain=c("Gt-19d1", "Gt-8d", "Gt-23d", "Gt-LH10", "Gt-4e",
                         "Gt-CB1", "Gt-3aA1", "Gh-2C17", "Gh-1B17"),
                file1=c("Gt17LH19d1", "Gt17LH48d", "Gt23d", "Gt14LH10", "Gt4e",
                        "GtCB1", "PG3aA1", "NZ1292C17", "Gh1B17"),
                file2=c("Gt-19d1", "Gt-8d", "Gt-23d", "Gt14LH10", "Gt-4e",
                        "Gt-CB1", "Gt-3aA1", "Gh-2C17", "Gh-1B17"),
                genespace=c("Gt19d1", "Gt8d", "Gt23d", "Gt14LH10", "Gt4e",
                            "GtCB1", "Gt3aA1", "Gaeumannomyces_hyphopodioides_2C17", "Gh1B17"))

plot_riparian_edit <- function(gsParam,
                               genomeIDs = gsParam$genomeIDs,
                               refGenome = NULL,
                               highlightBed = NULL,
                               blk = NULL,
                               backgroundColor = add_alpha("grey20", .5),
                               pdfFile = NULL,
                               useRegions = TRUE,
                               labelTheseGenomes = gsParam$genomeIDs,
                               useOrder = TRUE,
                               gapProp = 0.005,
                               chrFill = "white",
                               scalePlotHeight = 1,
                               scalePlotWidth = 1,
                               minChrLen2plot = ifelse(useOrder, 100, 1e5),
                               braidAlpha = .5,
                               scaleBraidGap = 0,
                               reorderBySynteny = TRUE,
                               howSquare = 0,
                               customRefChrOrder = NULL,
                               palette = gs_colors,
                               chrLabFontSize = 5,
                               chrExpand = 0.5,
                               chrBorderCol = chrFill,
                               syntenyWeight = .5,
                               chrBorderLwd = 0.2,
                               inversionColor = NULL,
                               invertTheseChrs = NULL,
                               forceRecalcBlocks = TRUE,
                               chrLabFun = function(x)
                                 gsub("^0", "",
                                      gsub("chr|scaf|chromosome|scaffold|^lg|_", "", tolower(x))),
                               xlabel = sprintf(
                                 "Chromosomes scaled by %s",
                                 ifelse(useOrder, "gene rank order", "physical position")),
                               scaleGapSize = 0.25,
                               addThemes = NULL,
                               verbose = FALSE,
                               refChrOrdFun = function(x)
                                 frank(list(as.numeric(gsub('\\D+','', x)), x), ties.method = "random")){
  
  color <- blkID <- pass <- genome <- chr <- start <- end <- index <- hasAll <-
    genome1 <- genome2 <- NULL
  # -- there are two methods of plotting.
  
  # 1. Default phase the blocks by a reference genome chromosomes
  # 2. A background of one color and then highlighted colors overlapping it.
  
  if(!"synteny" %in% names(gsParam))
    gsParam <- set_syntenyParams(gsParam)
  
  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
  
  if(is.null(highlightBed)){
    
    hapGs <- names(gsParam$ploidy)[gsParam$ploidy == 1]
    if(length(hapGs) == 0)
      stop("riparian plots must have a haploid reference as an anchor. Either re-run GENESPACE with a haploid outgroup or use custom parameters in `riparian_engine`")
    if(is.null(refGenome))
      refGenome <- hapGs[1]
    if(!refGenome %in% hapGs)
      stop(sprintf("specified reference genome %s is not haploid, but haploid genomes %s are in this GENESPACE run. Either use one of the haploid references or use custom parameters in `riparian_engine`",
                   refGenome, paste(hapGs, collapse = ", ")))
    
    bf <- file.path(gsParam$paths$riparian, sprintf("%s_phasedBlks.csv", refGenome))
    if(file.exists(bf) && !forceRecalcBlocks){
      blksTp <- fread(bf)
    }else{
      blksTp <- phase_blks(
        gsParam = gsParam,
        refGenome = refGenome,
        useRegions = useRegions,
        blkSize = gsParam$params$blkSize,
        synBuff = gsParam$params$synBuff,
        refChr = NULL,
        refStartBp = NULL,
        refEndBp = NULL)
    }
    
    if(is.null(pdfFile)){
      pltDat <- riparian_engine_edit(
        blk = blksTp, bed = bed, refGenome = refGenome,
        genomeIDs = genomeIDs, labelTheseGenomes = labelTheseGenomes,
        useOrder = useOrder, gapProp = gapProp, chrFill = chrFill,
        minChrLen2plot = minChrLen2plot, chrLabFontSize = chrLabFontSize,
        braidAlpha = braidAlpha, scaleBraidGap = scaleBraidGap, verbose = verbose,
        reorderBySynteny = reorderBySynteny, howSquare = howSquare,
        customRefChrOrder = customRefChrOrder, addThemes = addThemes,
        invertTheseChrs = invertTheseChrs, chrLabFun = chrLabFun, xlabel = xlabel,
        chrExpand = chrExpand, scaleGapSize = scaleGapSize, palette = palette,
        chrBorderCol = chrBorderCol, chrBorderLwd = chrBorderLwd,
        syntenyWeight = syntenyWeight,
        inversionColor = inversionColor, refChrOrdFun = refChrOrdFun)
      print(pltDat$ggplotObj)
    }else{
      pdf(
        file = pdfFile,
        height = sqrt(length(genomeIDs)) * scalePlotHeight * 2,
        width = 8 * scalePlotWidth)
      pltDat <- riparian_engine_edit(
        blk = blksTp, bed = bed, refGenome = refGenome,
        genomeIDs = genomeIDs, labelTheseGenomes = labelTheseGenomes,
        useOrder = useOrder, gapProp = gapProp, chrFill = chrFill,
        minChrLen2plot = minChrLen2plot, chrLabFontSize = chrLabFontSize,
        braidAlpha = braidAlpha, scaleBraidGap = scaleBraidGap, verbose = verbose,
        reorderBySynteny = reorderBySynteny, howSquare = howSquare,
        customRefChrOrder = customRefChrOrder, addThemes = addThemes,
        invertTheseChrs = invertTheseChrs, chrLabFun = chrLabFun, xlabel = xlabel,
        chrExpand = chrExpand, scaleGapSize = scaleGapSize, palette = palette,
        chrBorderCol = chrBorderCol, chrBorderLwd = chrBorderLwd,
        syntenyWeight = syntenyWeight,
        inversionColor = inversionColor, refChrOrdFun = refChrOrdFun)
      print(pltDat$ggplotObj)
      dev.off()
    }
  }else{
    highlightBed <- data.table(highlightBed)
    if(nrow(highlightBed) == 0)
      stop("highlightBed must be a data.table or data.frame with >= 1 rows\n")
    if(!all(c("chr", "genome") %in% names(highlightBed)))
      stop("highlightBed must be a data.table or data.frame with columns named chr and genome\n")
    u <- unique(paste(bed$genome, bed$chr))
    highlightBed[,pass := paste(genome, chr) %in% u]
    if(!all(highlightBed$pass))
      stop("the following genome/chr combinations are not in the run:", subset(bed, !pass)[,1:2])
    highlightBed[,pass := NULL]
    if(!"start" %in% names(highlightBed))
      highlightBed[,start := 0]
    if(!"end" %in% names(highlightBed))
      highlightBed[,end := Inf]
    
    if(is.null(backgroundColor) || is.na(backgroundColor)){
      allBlks <- NULL
    }else{
      hapGs <- names(gsParam$ploidy)[gsParam$ploidy == 1]
      if(is.null(refGenome))
        refGenome <- hapGs[1]
      
      # -- get background all blocks
      allBlks <- nophase_blks(
        gsParam = gsParam,
        useRegions = useRegions,
        blkSize = gsParam$params$blkSize)
      allBlks[,color := backgroundColor]
      allBlks[,blkID := sprintf("all_%s", blkID)]
    }
    
    # -- get blocks in the intervals of interest
    
    if(!"color" %in%  names(highlightBed))
      highlightBed[,color := gs_colors(.N)]
    blksBed <- rbindlist(lapply(1:nrow(highlightBed), function(i){
      suppressWarnings(blksTp <- with(highlightBed, phase_blks(
        gsParam = gsParam,
        refGenome = genome[i],
        useRegions = useRegions,
        blkSize = gsParam$params$blkSize,
        synBuff = gsParam$params$synBuff,
        refChr = chr[i],
        refStartBp = start[i],
        refEndBp = end[i])))
      blksTp[,color := highlightBed$color[i]]
      blksTp[,index := i]
      return(blksTp)
    }), fill = T)
    
    if(!"genome1" %in% colnames(blksBed))
      stop("none of the lines in highlightbed contain syntenic anchor genes\n")
    
    bads <- subset(blksBed, !complete.cases(blksBed))
    if(nrow(bads) > 0)
      warning(sprintf(
        "lines %s in highlight bed do not have any syntenic genes; these will be ignored",
        paste(unique(bads$index), collapse = ",")))
    goods <- subset(blksBed, complete.cases(blksBed))
    goods[,hasAll := all(genomeIDs %in% c(genome1, genome2)), by = "index"]
    bads <- subset(goods, !hasAll)
    if(nrow(bads) > 0)
      warning(sprintf(
        "lines %s in highlight bed do not connect all genomes, these may look strange in the plot",
        paste(unique(bads$index), collapse = ",")))
    blksTp <- rbind(allBlks, goods, fill = T)
    
    if(is.null(pdfFile)){
      pltDat <- riparian_engine_edit(
        blk = blksTp, bed = bed, refGenome = refGenome,
        theseBlocksFirst = unique(allBlks$blkID), verbose = verbose,
        genomeIDs = genomeIDs, labelTheseGenomes = labelTheseGenomes,
        useOrder = useOrder, gapProp = gapProp, chrFill = chrFill,
        minChrLen2plot = minChrLen2plot, chrLabFontSize = chrLabFontSize,
        braidAlpha = braidAlpha, scaleBraidGap = scaleBraidGap,
        reorderBySynteny = reorderBySynteny, howSquare = howSquare,
        customRefChrOrder = customRefChrOrder, addThemes = addThemes,
        invertTheseChrs = invertTheseChrs, chrLabFun = chrLabFun, xlabel = xlabel,
        scaleGapSize = scaleGapSize, refChrOrdFun = refChrOrdFun,
        chrBorderCol = chrBorderCol, chrBorderLwd = chrBorderLwd,
        inversionColor = inversionColor, syntenyWeight = syntenyWeight,
        chrExpand = chrExpand, palette = colorRampPalette(backgroundColor))
      print(pltDat$ggplotObj)
    }else{
      pdf(
        file = pdfFile,
        height = (length(genomeIDs) / 2) * scalePlotHeight,
        width = 8 * scalePlotWidth)
      pltDat <- riparian_engine_edit(
        blk = blksTp, bed = bed, refGenome = refGenome,
        theseBlocksFirst = unique(allBlks$blkID), verbose = verbose,
        genomeIDs = genomeIDs, labelTheseGenomes = labelTheseGenomes,
        useOrder = useOrder, gapProp = gapProp, chrFill = chrFill,
        minChrLen2plot = minChrLen2plot, chrLabFontSize = chrLabFontSize,
        braidAlpha = braidAlpha, scaleBraidGap = scaleBraidGap,
        reorderBySynteny = reorderBySynteny, howSquare = howSquare,
        customRefChrOrder = customRefChrOrder, addThemes = addThemes,
        invertTheseChrs = invertTheseChrs, chrLabFun = chrLabFun, xlabel = xlabel,
        scaleGapSize = scaleGapSize, refChrOrdFun = refChrOrdFun,
        chrBorderCol = chrBorderCol, chrBorderLwd = chrBorderLwd,
        inversionColor = inversionColor, syntenyWeight = syntenyWeight,
        chrExpand = chrExpand, palette = colorRampPalette(backgroundColor))
      print(pltDat$ggplotObj)
      dev.off()
    }
  }
  return(list(blks = blksTp, plotData = pltDat))
}


riparian_engine_edit <- function(blk,
                                 bed,
                                 refGenome = NULL,
                                 genomeIDs = NULL,
                                 theseBlocksFirst = NULL,
                                 labelTheseGenomes = NULL,
                                 useOrder = TRUE,
                                 gapProp = 0.005,
                                 chrFill = "white",
                                 chrLabFontSize = 5,
                                 minChrLen2plot = ifelse(useOrder, 100, 1e5),
                                 braidAlpha = .5,
                                 scaleBraidGap = 0,
                                 reorderBySynteny = TRUE,
                                 howSquare = 0,
                                 chrExpand = .5,
                                 chrBorderCol = NA,
                                 chrBorderLwd = 0,
                                 customRefChrOrder = NULL,
                                 palette = gs_colors,
                                 inversionColor = NULL,
                                 invertTheseChrs = NULL,
                                 syntenyWeight = .5,
                                 chrLabFun = function(x)
                                   gsub("^0", "",
                                        gsub("chr|scaf|chromosome|scaffold|^lg|_", "", tolower(x))),
                                 xlabel = sprintf(
                                   "Chromosomes scaled by %s",
                                   ifelse(useOrder, "gene rank order", "physical position")),
                                 scaleGapSize = 0.25,
                                 addThemes = NULL,
                                 verbose = FALSE,
                                 refChrOrdFun = function(x)
                                   frank(list(as.numeric(gsub('\\D+','', x)), x), ties.method = "random")){
  
  ##############################################################################
  subset_toChainedGenomes <- function(genomeIDs, blk){
    genome1 <- genome2 <- NULL
    genomeOrd <- data.table(
      genome1 = genomeIDs[-length(genomeIDs)], genome2 = genomeIDs[-1])
    genomeOrd[,`:=`(y1 = match(genome1, genomeIDs),
                    y2 = match(genome2, genomeIDs))]
    u <- with(genomeOrd, paste(genome1, genome2))
    if(!all(u %in% paste(blk$genome1, blk$genome2)))
      warning("problem with the block coordinates, some chained combinations of genomeIDs are not in the blocks\n")
    blk <- merge(genomeOrd, blk, by = c("genome1", "genome2"))
    return(blk)
  }
  
  ##############################################################################
  get_chrLens <- function(blk, minChrLen2plot, useOrder){
    genome1 <- genome2 <- chr1 <- chr2 <- startBp1 <- startBp2 <- endBp1 <-
      endBp2 <- startOrd1 <- startOrd2 <- endOrd1 <- endOrd2 <- end <- length <-
      genome <- chr <- NULL
    if(useOrder){
      chrLengths <- with(blk, data.table(
        genome = c(genome1, genome2, genome1, genome2),
        chr = c(chr1, chr2, chr1, chr2),
        end = c(startOrd1, startOrd2, endOrd1, endOrd2)))
    }else{
      chrLengths <- with(blk, data.table(
        genome = c(genome1, genome2, genome1, genome2),
        chr = c(chr1, chr2, chr1, chr2),
        end = c(startBp1, startBp2, endBp1, endBp2)))
    }
    
    chrLengths <- chrLengths[,list(length = max(end)), by = c("genome", "chr")]
    u <- with(blk, unique(paste(c(genome1, genome2), c(chr1, chr2))))
    if(!all(u %in% paste(chrLengths$genome, chrLengths$chr))){
      warning("could not parse chrLengths, some genome/chr combinations in block coordinates are not in chrLengths\n")
      chrLengths <- NULL
    }
    
    out <- subset(chrLengths, length >= minChrLen2plot)
    return(out)
  }
  
  ##############################################################################
  make_scaleBarData <- function(chrLengths, ylabBuff, useOrder){
    length <- genome <- s <- chrstart <- genome <- bot <- top <- mid <- NULL
    s <- min(with(chrLengths, tapply(length, genome, sum))) * .3
    if(useOrder){
      s <- ifelse(s > 1e5, round(s, -5),
                  ifelse(s > 1e4, round(s, -4),
                         ifelse(s > 1e3, round(s, -3),
                                ifelse(s > 100, round(s, -2), round(s, -1)))))
    }else{
      s <- ifelse(s > 1e9, round(s, -9),
                  ifelse(s > 1e8, round(s, -8),
                         ifelse(s > 1e7, round(s, -7),
                                ifelse(s > 1e6, round(s, -6), round(s, -5)))))
    }
    
    top <- max(chrLengths$y2) + (ylabBuff * 3)
    bot <- max(top) - ylabBuff
    left <- max(with(chrLengths, tapply(chrstart, genome, min)))
    right <- left + s
    mid <- (bot + top)/2
    scbar <- data.table(
      line = c("left", "right","mid"),
      x = c(left, right, left),
      xend = c(left, right, right),
      y = c(bot, bot, mid),
      yend = c(top, top, mid))
    
    sclab <- sprintf(
      "%s %s", ifelse(useOrder, s, s/1e6), ifelse(useOrder, "genes", "Mbp"))
    
    return(list(label = sclab, segments = scbar))
  }
  
  ##############################################################################
  pull_synChrOrd <- function(refGenome, bed, clens, syntenyWeight){
    
    ct_clust <- function(chr, ord, nclus){
      if(any(duplicated(chr)))
        stop("can't have duplicated chr values")
      
      if(length(ord) == 1){
        names(ord) <- chr[1]
        return(ord)
      }else{
        if(length(ord) == 0){
          return(1)
        }else{
          y <- data.matrix(ord)
          rownames(y) <- chr
          di <- dist(y)
          hc <- hclust(di, method = "ward.D2")
          nclus <- min(c(length(ord), nclus))
          ct <- cutree(hc, k = nclus)
          return(ct[chr])
        }
      }
    }
    
    ordByFun <- genome <- noAnchor <- refchr <- reford <- ord <- median <-
      plotOrd <- meanPos <- isArrayRep <- refChrOrder <- pSyn <- pName <-
      synPos <- meanOrd <- NULL
    
    tb <- subset(bed, !noAnchor & isArrayRep)[,c("genome", "chr", "ord", "og")]
    tc <- clens[,c("genome", "chr", "ordByFun")]
    bedi <- merge(
      subset(tb, !duplicated(tb)),
      subset(tc, !duplicated(tc)),
      by = c("genome", "chr"), allow.cartesian = T)
    
    bedr <- with(subset(bedi, genome == refGenome), data.table(
      refgenome = refGenome, refchr = chr, refChrOrder = ordByFun,
      reford = ord, og = og))
    beda <- subset(bedi, genome != refGenome)
    
    bedm <- merge(
      subset(bedr, !duplicated(bedr)),
      subset(beda, !duplicated(beda)),
      by = "og", allow.cartesian = T)
    setkey(bedm, refChrOrder, reford, genome, ord)
    bedm[,reford := (1:.N)/.N, by = "genome"]
    bedm[,`:=`(pName = ordByFun / max(ordByFun),
               pSyn = reford/max(reford)), by = "genome"]
    
    chro <- bedm[,list(
      meanPos = median(pSyn, na.rm = T),
      meanOrd = median(pName, na.rm = T)),
      by = c("genome", "chr")]
    
    if(syntenyWeight <= 0){
      nclus <- 1
    }else{
      if(syntenyWeight >= 1){
        nclus <- Inf
      }else{
        nclus <- uniqueN(bedm$refchr)
      }
    }
    
    chro[,`:=`(synClus = ct_clust(ord = meanPos, chr = chr, nclus = nclus)),
         by = "genome"]
    chro[,`:=`(synPos = mean(meanPos, na.rm = T)), by = c("genome", "synClus")]
    setkey(chro, genome, synPos, meanOrd, chr)
    
    outa <- merge(clens, chro, by = c("genome", "chr"))
    outr <- subset(clens, genome == refGenome)
    outr[,plotOrd := ordByFun]
    
    setkey(outa, synPos, ordByFun, genome)
    outa[,plotOrd := 1:.N, by = "genome"]
    
    outa <- outa[,colnames(outr), with = F]
    
    return(rbind(outr, outa))
  }
  
  ##############################################################################
  reorder_inChrs <- function(blk){
    genome1 <- genome2 <- chr1 <- chr2 <- startOrd1 <- startOrd2 <- start <-
      genome <- chr <- endOrd1 <- endOrd2 <- chr1st <- chr2st <- mpv <- NULL
    minPos <- with(blk, data.table(
      genome = c(genome1, genome2),
      chr = c(chr1, chr2),
      start = c(startOrd1, startOrd2)))
    minPos <- minPos[,list(chrStart = min(start)), by = c("genome", "chr")]
    mpv <- minPos$chrStart; names(mpv) <- with(minPos, paste(genome, chr))
    blk[,`:=`(chr1st = mpv[paste(genome1, chr1)],
              chr2st = mpv[paste(genome2, chr2)])]
    blk[,`:=`(startOrd1 = (startOrd1 - chr1st) + 1,
              startOrd2 = (startOrd2 - chr2st) + 1,
              endOrd1 = (endOrd1 - chr1st) + 1,
              endOrd2 = (endOrd2 - chr2st) + 1)]
    return(blk)
  }
  
  ##############################################################################
  calc_gapSize <- function(chrLens, scaleGapSize, gapProp, rg){
    maxGapSize <- genomeSize <- totLen <- ngaps <- totLen <- genome <-
      diffMax <- gapSize <- NULL
    maxGapSize <- round(gapProp * sum(chrLens$length[chrLens$genome == rg]))
    genomeSize <- chrLens[,list(totLen = sum(as.numeric(length)), ngaps = .N),
                          by = "genome"]
    genomeSize[,diffMax := (max(totLen) - totLen)/ngaps]
    genomeSize[,gapSize := maxGapSize + (diffMax * scaleGapSize)]
    gapSize <- genomeSize$gapSize; names(gapSize) <- genomeSize$genome
    return(gapSize)
  }
  
  
  blkID <- refChr <- genome1 <- genome2 <- genome <- ordByFun <- chr <-
    plotOrd <- gapSize <- chrstart <- med <- ht <- y1 <- y2 <- chrLab <- u1 <-
    xs1 <- cstart <- xe1 <- cend <- difstart <- difend <- u2 <- xs2 <- xe2 <-
    x <- y <- color <- x1 <- x2 <- xend <- yend <- line <- NULL
  
  ##############################################################################
  ##############################################################################
  refG <- refGenome
  gids <- genomeIDs
  
  blk <- data.table(blk)
  blksFileOrd <- unique(blk$blkID)
  refGenome <- NULL
  genomeIDs <- NULL
  
  if(is.null(gids))
    gids <- unique(unlist(blk[,c("genome1", "genome2")]))
  fullblk <- data.table(blk)
  if(is.null(refG))
    refG <- gids[1]
  if(is.null(labelTheseGenomes))
    labelTheseGenomes <- gids
  
  if(!is.null(invertTheseChrs)){
    if(!is.data.frame(invertTheseChrs))
      stop("invertTheseChrs given, but this needs to be a data.frame\n")
    if(!all(c("genome", "chr") %in% colnames(invertTheseChrs)))
      stop("invertTheseChrs data.frame given, but does not have col names 'chr' and 'genome' \n")
    invchr <- with(invertTheseChrs, paste(genome, chr))
  }else{
    invchr <- NULL
  }
  
  if(!is.null(customRefChrOrder)){
    customRefChrOrder <- as.character(customRefChrOrder)
    if(any(is.null(customRefChrOrder)) || any(is.na(customRefChrOrder)))
      stop("customRefChrOrder is given but can't be coerced to character\n")
  }
  
  ##############################################################################
  # 1 Get chromosome positions nailed down
  # -- 1.1 if useorder, re-scale order so that each chr starts at 1
  if("refChr" %in% colnames(blk))
    blk[,blkID := paste(blkID, refChr)]
  blk <- subset(blk, genome1 %in% gids & genome2 %in% gids)
  if(useOrder)
    blk <- reorder_inChrs(blk)
  
  # -- 1.2 calculate chromosome lengths
  clens <- get_chrLens(
    blk = blk,
    useOrder = useOrder,
    minChrLen2plot = minChrLen2plot)
  
  # -- 1.3 get chromosome order by chromosome name
  if(!is.null(customRefChrOrder)){
    rlen <- subset(clens, genome == refG)
    if(any(!rlen$chr %in% customRefChrOrder))
      stop("customRefChrOrder but not all refchrs are in this list\n")
    rlen[,ordByFun := as.integer(factor(chr, levels = customRefChrOrder))]
    nlen <- subset(clens, genome != refG)
    nlen[,ordByFun := refChrOrdFun(chr), by = "genome"]
    clens <- rbind(rlen, nlen)
  }else{
    clens[,ordByFun := refChrOrdFun(chr), by = "genome"]
  }
  setkey(clens, genome, ordByFun)
  
  # -- 1.4 [optional] replace non-reference genome chr order to max synteny
  if(reorderBySynteny){
    # refGenome <<- refG
    # bed <<- bed
    # clens <<- clens
    # syntenyWeight <<- syntenyWeight
    clens <- pull_synChrOrd(
      refGenome = refG, bed = bed, clens = clens, syntenyWeight = syntenyWeight)
  }else{
    clens[,plotOrd := ordByFun]
  }
  clens[,ordByFun := NULL]
  
  # -- 1.5 get genome-specific gap sizes
  setkey(clens, genome, plotOrd)
  gapSizes <- calc_gapSize(
    chrLens = clens, scaleGapSize = scaleGapSize,
    gapProp = gapProp, rg = refG)
  
  # -- 1.6 add gap to the chrlengths
  clens[,gapSize := gapSizes[genome]]
  
  # -- 1.7 add cumulative position to the chromosomes
  clens[,chrstart := cumsum(c(1, (length[-.N]+1)) + gapSize)-gapSize,
        by = "genome"]
  
  # -- 1.8 get median position for each genome
  clens[,med := max(chrstart + length)/2,
        by = "genome"]
  
  # -- 1.9 add start and end positions of each chromosome in linear coords
  clens[,`:=`(x1 = chrstart - med, x2 = (chrstart + length) - med)]
  
  # -- 1.10.5 get plotting info so that we scale the chrs correctly
  chrLabFontHeight <- chrLabFontSize/72 # inches (real life)
  
  yrange <- length(gids) + .5
  ds <- dev.size(units = "in")
  plotWidth <- ds[1]
  plotHeight <- ds[2]
  
  nwords <- plotHeight / chrLabFontHeight
  wordHt <-  yrange / nwords
  chrPolyHeight <- wordHt + ((chrExpand * wordHt) / 2)  # inches (plotting
  
  # -- 1.10 add thickness of the chr polygons
  clens[,ht := ifelse(genome %in% labelTheseGenomes, chrPolyHeight/2, chrPolyHeight/4)]
  clens[,`:=`(y1 = match(genome, gids),
              y2 = match(genome, gids))]
  
  # -- 1.11 pull vector of chromosome start positions
  chrstartv <- clens$x1
  chrytop <- clens$y1 + ((chrPolyHeight / 2) * scaleBraidGap)
  chrybot <- clens$y2 - ((chrPolyHeight / 2) * scaleBraidGap)
  
  names(chrstartv) <- names(chrybot) <-
    names(chrytop) <- with(clens, paste(genome, chr))
  
  clens[,`:=`(y1 = y1 + ht,
              y2 = y2 - ht,
              ht = NULL)]
  
  # -- 1.12 make the polygons
  chrPolygons <- rbindlist(lapply(1:nrow(clens), function(i){
    z <- clens[i,]
    out <- data.table(z[,c("genome","chr")], with(z, round_rect(
      xleft = x1, xright = x2, ybottom = y2,
      ytop = y1, yrange = range(c(clens$y1, clens$y2)),
      xrange = range(c(clens$x1, clens$x2)),
      plotWidth = plotWidth + (plotWidth * howSquare),
      plotHeight = plotHeight)))
    return(out)
  }))
  
  
  # -- 1.13 add chromosome labels
  clens[,chrLab := chrLabFun(chr)]
  if(!is.null(invchr)){
    wh <- which(paste(clens$genome, clens$chr) %in% invchr)
    clens$chrLab[wh] <- sprintf("%s*", clens$chrLab[wh])
  }
  
  ##############################################################################
  # 2 Get block positions nailed down
  
  # -- 2.1 subset blocks to daisy chains
  dsychn <- subset_toChainedGenomes(genomeIDs = gids, blk = blk)
  
  # -- 2.2 add color column
  if("color" %in% colnames(dsychn)){
    if(any(!are_colors(dsychn$color)))
      stop("colors are given in the blocks but can't be coerced to R colors\n")
  }else{
    ulevs <- subset(clens, genome == dsychn$refGenome[1])[,c("plotOrd", "chr")]
    setkey(ulevs, plotOrd)
    ulevs[,`:=`(color = palette(.N), refChr = chr, plotOrd = NULL, chr = NULL)]
    dsychn <- merge(dsychn, ulevs, by = "refChr")
  }
  
  if(!is.null(inversionColor)){
    if(are_colors(inversionColor[1]))
      dsychn$color[dsychn$orient == "-"] <- inversionColor[1]
  }
  
  # -- 2.4 simplify and convert to projections depending on useOrder
  if(useOrder){
    tp <- with(dsychn, data.table(
      xs1 = startOrd1 + chrstartv[paste(genome1, chr1)],
      xs2 = startOrd2 + chrstartv[paste(genome2, chr2)],
      xe1 = endOrd1 + chrstartv[paste(genome1, chr1)],
      xe2 = endOrd2 + chrstartv[paste(genome2, chr2)],
      y1 = chrytop[paste(genome1, chr1)],
      y2 = chrybot[paste(genome2, chr2)],
      u1 = paste(genome1, chr1),
      u2 = paste(genome2, chr2),
      color = color, blkID = blkID))
  }else{
    tp <- with(dsychn, data.table(
      xs1 = startBp1 + chrstartv[paste(genome1, chr1)],
      xs2 = startBp2 + chrstartv[paste(genome2, chr2)],
      xe1 = endBp1 + chrstartv[paste(genome1, chr1)],
      xe2 = endBp2 + chrstartv[paste(genome2, chr2)],
      y1 = chrytop[paste(genome1, chr1)],
      y2 = chrybot[paste(genome2, chr2)],
      u1 = paste(genome1, chr1),
      u2 = paste(genome2, chr2),
      color = color, blkID = blkID))
  }
  
  # -- 2.5 if necessary invert some chromosomes
  if(!is.null(invchr)){
    tmp <- subset(clens, paste(genome, chr) %in% invchr)
    chrStarts <- tmp$x1; names(chrStarts) <- with(tmp, paste(genome, chr))
    chrEnds <- tmp$x2; names(chrEnds) <- with(tmp, paste(genome, chr))
    
    tp1 <- subset(tp, u1 %in% invchr)
    if(nrow(tp1) > 0){
      tp0 <- subset(tp, !u1 %in% invchr)
      tp1[,`:=`(cstart = chrStarts[u1], cend = chrEnds[u1])]
      tp1[,`:=`(difstart = xs1 - cstart, difend = xe1 - cstart)]
      tp1[,`:=`(xs1 = cend - difstart, xe1 = cend - difend)]
      tp1[,`:=`(cstart = NULL, cend = NULL, difstart = NULL, difend = NULL)]
      tp <- rbind(tp0, tp1)
    }
    
    tp2 <- subset(tp, u2 %in% invchr)
    if(nrow(tp2) > 0){
      tp0 <- subset(tp, !u2 %in% invchr)
      tp2[,`:=`(cstart = chrStarts[u2], cend = chrEnds[u2])]
      tp2[,`:=`(difstart = xs2 - cstart, difend = xe2 - cstart)]
      tp2[,`:=`(xs2 = cend - difstart, xe2 = cend - difend)]
      tp2[,`:=`(cstart = NULL, cend = NULL, difstart = NULL, difend = NULL)]
      tp <- rbind(tp0, tp2)
    }
  }
  
  # ADD TELOMERES
  
  for (i in 1:length(strains$strain)) {
    
    telomeres <- read.csv(
      paste0("/ei/.project-scratch/d/d2c0bfb1-c37a-4211-9493-86b15d4e773e/015_tapestry/",
             strains$file1[i], "/", "contig_details.tsv"),
      sep="\t"
    )
    
    telomeres$genome <- strains$strain[i]
    
    telomeres$Contig <- paste0(strains$genespace[i], "_EI_v1.1_",
                               sub("l", "", telomeres$Contig))
    
    telomeres$StartTelomeres <- ifelse(telomeres$StartTelomeres >= 5, TRUE, FALSE)
    telomeres$EndTelomeres <- ifelse(telomeres$EndTelomeres >= 5, TRUE, FALSE)
    
    assign(paste0("telomeres.", strains$strain[i]), telomeres)
    
  }
  
  all.telomeres <- do.call("rbind", mget(ls(pattern="telomeres.")))
  
  all.telomeres <- all.telomeres[which(all.telomeres$StartTelomeres == TRUE |
                                         all.telomeres$EndTelomeres == TRUE),]
  
  all.telomeres <- all.telomeres[all.telomeres$Contig %in% unique(chrPolygons$chr),]
  
  all.telomeres[,c("start.x", "start.y", "start.yend",
                   "end.x", "end.y", "end.yend")] <- NA
  
  for (i in 1:length(all.telomeres$Contig)) {
    
    if (all.telomeres$Contig[i] %in% invertTheseChrs$chr) {
      
      left <- all.telomeres$StartTelomeres[i]
      right <- all.telomeres$EndTelomeres[i]
      
      all.telomeres$EndTelomeres[i] <- left
      all.telomeres$StartTelomeres[i] <- right
      
    }
    
    tmp <- chrPolygons[chrPolygons$chr == all.telomeres$Contig[i],]
    
    if (all.telomeres$StartTelomeres[i]) {
      
      all.telomeres$start.x[i] <- min(tmp$x)
      all.telomeres$start.y[i] <- min(tmp$y)
      all.telomeres$start.yend[i] <- max(tmp$y)
      
    }
    
    if (all.telomeres$EndTelomeres[i]) {
      
      all.telomeres$end.x[i] <- max(tmp$x)
      all.telomeres$end.y[i] <- min(tmp$y)
      all.telomeres$end.yend[i] <- max(tmp$y)
      
    }
    
  }
  
  # -- 2.6 make the braid polygons
  braidPolygons <- rbindlist(lapply(1:nrow(tp), function(i){
    x <- with(tp[i, ], calc_curvePolygon(
      start1 = xs1, end1 = xe1, start2 = xs2,
      end2 = xe2, y1 = y1, y2 = y2))
    x[,`:=`(blkID = tp$blkID[i], color = tp$color[i])]
    return(x)
  }))
  
  # -- 6. Make the scalebar
  sbData <- make_scaleBarData(
    chrLengths = clens, ylabBuff = chrPolyHeight, useOrder = useOrder)
  
  glab <- subset(clens, !duplicated(genome))
  setkey(glab, y1)
  
  if(is.null(theseBlocksFirst)){
    p <- ggplot()+
      
      # -- all blocks
      geom_polygon(
        data = braidPolygons,
        aes(x = x, y = y, group = blkID, fill = color),
        alpha = braidAlpha, lwd = NA)
  }else{
    p <- ggplot()+
      
      # -- bg blocks
      geom_polygon(
        data = subset(braidPolygons, blkID %in% theseBlocksFirst),
        aes(x = x, y = y, group = blkID, fill = color),
        alpha = braidAlpha, lwd = NA)+
      
      # -- roi blocks
      geom_polygon(
        data = subset(braidPolygons, !blkID %in% theseBlocksFirst),
        aes(x = x, y = y, group = blkID, fill = color),
        alpha = braidAlpha, lwd = NA)
  }
  p <- p +
    
    # -- braid colors
    scale_fill_identity(guide = "none")+
    
    # -- chromosome label backgrounds
    
    geom_polygon(
      data = chrPolygons,
      aes(x = x, y = y, group = paste(genome, chr)),
      fill = chrFill, lwd = chrBorderLwd, col = chrBorderCol)+
    
    # -- chr label text
    geom_label(
      data = subset(clens, genome %in% labelTheseGenomes),
      aes(x = (x1+x2)/2, y = (y1+y2)/2, label = chrLab),
      col = "black", size = chrLabFontSize*0.36, hjust = .5,
      label.padding = unit(0.1, "lines"), label.size = 0,
      alpha = 0.5, nudge_y = 0.3)+
    
    # -- add start telomeres
    geom_segment(
      data = all.telomeres,
      aes(x = start.x,
          xend = start.x,
          y = start.y-0.03,
          yend = start.yend+0.03),
      colour = "black",
      linewidth = 0.5) +
    # -- add end telomeres
    geom_segment(
      data = all.telomeres,
      aes(x = end.x,
          xend = end.x,
          y = end.y-0.03,
          yend = end.yend+0.03),
      colour = "black",
      linewidth = 0.5)+
    
    # -- scale bar
    geom_segment(
      data = sbData$segments,
      aes(x = x, xend = xend, y = y, yend = yend),
      lwd = chrBorderLwd, col = chrBorderCol)+
    geom_text(
      data = subset(sbData$segments, line == "mid"),
      aes(x = (x + xend)/2, y = y),
      label = sbData$label,
      col = chrBorderCol, size = chrLabFontSize*0.36, hjust = .5, vjust = 1.5)+
    
    # -- genome labels
    scale_y_continuous(breaks = (glab$y1+glab$y2)/2, labels = glab$genome, expand = c(0.01, 0.01), name = NULL)+
    scale_x_continuous(expand = c(0.005, 0.005), name = xlabel)+
    
    # -- themes
    theme(panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank())
  
  if(!is.null(addThemes))
    p <- p + addThemes
  
  # chrPolygons <<- chrPolygons
  # p <<- p
  
  #Add tree
  p <- p %>% insert_left(gg.tree, width=0.3)
  
  return(list(ggplotObj = p,
              sourceData = list(blocks = dsychn, chromosomes = clens)))
}

#GENESPACE run working directory
wd <- "/ei/.project-scratch/d/d2c0bfb1-c37a-4211-9493-86b15d4e773e/genespace_gaeumannomyces_230510"

#Dataframe of fragments to invert to optimise visualisation
inversions <-
  data.frame(genome=c("Gt8d",
                      "Gt8d",
                      "Gt8d",
                      "Gt8d",
                      
                      "Gt4e",
                      "Gt4e",
                      "Gt4e",
                      
                      "Gt23d",
                      
                      "Gt19d1",
                      "Gt19d1",
                      
                      "Gt14LH10",
                      "Gt14LH10",
                      "Gt14LH10",
                      
                      "Gt3aA1",
                      "Gt3aA1",
                      "Gt3aA1",
                      
                      "GtCB1",
                      "GtCB1",
                      "GtCB1",
                      "GtCB1",
                      "GtCB1",
                      "GtCB1",
                      
                      "Gh2C17",
                      "Gh2C17",
                      "Gh2C17",
                      
                      "Gh1B17",
                      "Gh1B17",
                      "Gh1B17",
                      "Gh1B17",
                      "Gh1B17",
                      "Gh1B17"),
             chr=c("Gt8d_EI_v1.1_ptg000005",
                   "Gt8d_EI_v1.1_ptg000001",
                   "Gt8d_EI_v1.1_ptg000002",
                   "Gt8d_EI_v1.1_ptg000003",
                   
                   "Gt4e_EI_v1.1_ptg000008",
                   "Gt4e_EI_v1.1_ptg000004",
                   "Gt4e_EI_v1.1_ptg000002",
                   
                   "Gt23d_EI_v1.1_ptg000006",
                   
                   "Gt19d1_EI_v1.1_ptg000003",
                   "Gt19d1_EI_v1.1_ptg000006",
                   
                   "Gt14LH10_EI_v1.1_ptg000004",
                   "Gt14LH10_EI_v1.1_ptg000005",
                   "Gt14LH10_EI_v1.1_ptg000010",
                   
                   "Gt3aA1_EI_v1.1_ptg000001",
                   "Gt3aA1_EI_v1.1_ptg000006",
                   "Gt3aA1_EI_v1.1_ptg000007",
                   
                   "GtCB1_EI_v1.1_ptg000002",
                   "GtCB1_EI_v1.1_ptg000005",
                   "GtCB1_EI_v1.1_ptg000007",
                   "GtCB1_EI_v1.1_ptg000001",
                   "GtCB1_EI_v1.1_ptg000008",
                   "GtCB1_EI_v1.1_ptg000004",
                   
                   "Gaeumannomyces_hyphopodioides_2C17_EI_v1.1_ptg000001",
                   "Gaeumannomyces_hyphopodioides_2C17_EI_v1.1_ptg000003",
                   "Gaeumannomyces_hyphopodioides_2C17_EI_v1.1_ptg000006",
                   
                   "Gh1B17_EI_v1.1_ptg000002",
                   "Gh1B17_EI_v1.1_ptg000003",
                   "Gh1B17_EI_v1.1_ptg000015",
                   "Gh1B17_EI_v1.1_ptg000005",
                   "Gh1B17_EI_v1.1_ptg000007",
                   "Gh1B17_EI_v1.1_ptg000018"
             ))

#Load GENESPACE run parameters
load(paste0(wd, "/results/gsParams.rda"), verbose=TRUE)

#Use 1 core
gsParam$params$nCores <- 1

#Customise plotting theme
theme_custom <- function(){
  
  theme_genespace() %+replace%
    
    theme(
      panel.background=element_blank(),
      axis.text.y=element_blank()
    )
  
}

#Read in pseudochromosome designation to rename fragments
pseudochromosomes <- read.csv("pseudochromosomes.tsv",
                              sep="\t",
                              header=FALSE)

#Plot GENESPACE riparian plot
plot_riparian_edit(chrBorderCol="dimgrey",
                   refGenome="Gt23d",
                   invertTheseChrs=inversions,
                   scalePlotHeight=0.5,
                   scalePlotWidth=7/8,
                   customRefChrOrder=c("Gt23d_EI_v1.1_ptg000003",
                                       "Gt23d_EI_v1.1_ptg000005",
                                       "Gt23d_EI_v1.1_ptg000006",
                                       "Gt23d_EI_v1.1_ptg000004",
                                       "Gt23d_EI_v1.1_ptg000001",
                                       "Gt23d_EI_v1.1_ptg000002",
                                       "Gt23d_EI_v1.1_ptg000007"),
                   reorderBySynteny=TRUE,
                   syntenyWeight=1,
                   palette=colorRampPalette(
                     c("#FFAABB", "#EE8866", "#EE8866", "#EEDD88",
                                "#BBCC33", "#44BB99", "#77AADD")),
                                braidAlpha=0.5,
                   chrLabFun=function(x)
                     sub("pseu_chr", "",
                         pseudochromosomes$V3[
                           match(x, sub(
                             "-", "", pseudochromosomes$V2
                           ))
                         ]),
                   howSquare=Inf,
                   gsParam=gsParam,
                   useRegions=FALSE,
                   pdfFile=paste0("gaeumannomyces_synteny_", Sys.Date(), ".pdf"),
                   addThemes=theme_custom())
