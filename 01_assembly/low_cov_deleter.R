#!/usr/bin/env Rscript

library(dplyr)
library(Biostrings)

readKatHist <- function(file) {
  text <- tail(scan(file = file, what = "character"), -19)
  x <- strtoi(text[c(T, F)])
  y <- strtoi(text[!c(T, F)])
  return(list(x = x, y = y))
}

build_persistent_homology <- function(x) {
  peaks <- list()
  idx2peak <- rep(-1, length(x))
  indicies <- order(x, decreasing = T)
  
  for(idx in indicies){
    lftdone <- idx > 1 && idx2peak[idx - 1] != -1
    rgtdone <- idx < length(x) && idx2peak[idx + 1] != -1
    il <- ifelse(lftdone, idx2peak[idx - 1], -1) 
    ir <- ifelse(rgtdone, idx2peak[idx + 1], -1)
  
    # New peak born.
    if(!lftdone && !rgtdone){
      peaks[[length(peaks) + 1]] <- list(born = idx, died = -1, left = idx, right = idx)
      idx2peak[idx] <- length(peaks)
    }

    # Directly merge to next peak left.
    if(lftdone && !rgtdone) {
      peaks[[il]]$right <- peaks[[il]]$right + 1
      idx2peak[idx] <- il
    }

    # Directly merge to next peak right.
    if(!lftdone && rgtdone) {
      peaks[[ir]]$left <- peaks[[ir]]$left - 1
      idx2peak[idx] <- ir
    }

    # Merge left and right peaks
    if(lftdone && rgtdone) {
      # Left was born earlier: merge right to left.
      if(x[peaks[[il]]$born] > x[peaks[[ir]]$born]) {
      peaks[[ir]]$died <- idx
      peaks[[il]]$right <- peaks[[ir]]$right
      idx2peak[peaks[[il]]$right] <- idx2peak[idx] <- il
      } else {
        peaks[[il]]$died <- idx
        peaks[[ir]]$left <- peaks[[il]]$left
        idx2peak[peaks[[ir]]$left] <- idx2peak[idx] <- ir
      }
    }
  }
  return(peaks)
}

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

find_valleys <- function(x, m = 3) { return(find_peaks(-x, m)) }

compute_hist_threshold <- function(histFile) {
  katHist <- readKatHist(histFile)
  histPersistentHomology <- build_persistent_homology(katHist$y)
  histPersistentHomology <- lapply(histPersistentHomology, function(x){
    x$persistence <- ifelse(x$died == -1, Inf, katHist$y[x$born] - katHist$y[x$died])
    return(x)
  })
  histPersistentHomology <- histPersistentHomology[order(unlist(lapply(histPersistentHomology, function(x) x$persistence)), decreasing = T)]
  histThreshold <- histPersistentHomology[[2]]$left
  return(list(hist = katHist, homology = histPersistentHomology, threshold = histThreshold))
}

compute_sect_threshold <- function(sectFile, histAnalysis) {
  katSect <- arrange(read.table(sectFile, header = T), desc(seq_length), rev = T)
  sectHist <- table(katSect$median)
  sectValleys <- find_valleys(sectHist)
  sectThresholds <- strtoi(names(sectValleys))
  return(list(sect = katSect, hist = sectHist, valleys = sectValleys, putativeThresholds = sectThresholds))
}

compute_final_threshold <- function(histAnalysis, sectAnalysis, chooser = max) {
  pThresh <- sectAnalysis$putativeThresholds
  histPermittedThresholds <- pThresh[pThresh <= histAnalysis$threshold]
  return(list(histPermittedThresholds = histPermittedThresholds, chosenThreshold = chooser(histPermittedThresholds)))
}

filter_assembly <- function(infile, sectAnalysis, histAnalysis) {
  assemblySequences  <- readDNAStringSet(infile, format = "fasta")
  keepSequences <- sectAnalysis$sect %>%
    arrange(seq_name) %>%
    filter(median > histAnalysis$threshold) %>%
    pull(seq_name)
  filteredSeqs <- assemblySequences[keepSequences]
  return(list(fullAssembly = assemblySequences, filteredAssembly = filteredSeqs))
}

analysisMain <- function(asmFile, katHistFile, katSectFile, outputFasta) {
  message(paste0("Filtering: ", asmFile))
  histAnalysis <- compute_hist_threshold(katHistFile)
  sectAnalysis <- compute_sect_threshold(katSectFile)
  message(str(sectAnalysis$hist))
  message(str(sectAnalysis$valleys))
  message(paste0("Sect: putative threshold: ", sectAnalysis$putativeThresholds, "\n"))
  message(paste0("Hist: threshold: ", histAnalysis$threshold, "\n"))
#  thresholdDescision <- compute_final_threshold(histAnalysis, sectAnalysis)
#  message(paste0("Chosen kmer-cov threshold: ", thresholdDescision$chosenThrreshold))
  assemblyAnalysis <- filter_assembly(asmFile, sectAnalysis, histAnalysis)
  message(paste0("Sequences in: ", length(assemblyAnalysis$fullAssembly), "\n"))
  message(paste0("Sequences out: ", length(assemblyAnalysis$filteredAssembly), "\n"))
  writeXStringSet(assemblyAnalysis$filteredAssembly, outputFasta, format="fasta")
  q(save="no")
}

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 4){
  stop("Bad number of arguments", call.=FALSE)
}

analysisMain(args[1], args[2], args[3], args[4])



