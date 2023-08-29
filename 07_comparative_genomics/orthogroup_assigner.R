library(tidyverse)
library(rvest)

## Read in data ##

#GENESPACE directory
dir <- "S:/genespace_gaeumannomyces_230510/orthofinder/Results_May10/"

#Read in orthogroups from OrthoFinder run
orthogroups <- read.csv(paste0(dir, "Phylogenetic_Hierarchical_Orthogroups/N0.tsv"),
                        row.names=1, sep="\t", check.names=FALSE)
orthogroups <- orthogroups[,-which(c("OG", "Gene Tree Parent Clade") %in% colnames(orthogroups))]

#Read in 'unassigned genes' - i.e. strain specific genes - and combine with orthogroups dataframe
unassigned <- read.csv(paste0(dir, "Orthogroups/Orthogroups_UnassignedGenes.tsv"),
                       row.names=1, sep="\t", check.names=FALSE)
orthogroups <- rbind(orthogroups, unassigned)

#Remove outgroup
orthogroups <- orthogroups[,-which(colnames(orthogroups) == "Magpoae")]
orthogroups <- orthogroups[apply(orthogroups, 1, function(x) any(x != "")), ]

#Fix column names
colnames(orthogroups) <- c("Gh-1B17", "Gh-2C17", "Gt14LH10", "Gt-19d1",
                           "Gt-23d", "Gt-3aA1", "Gt-4e", "Gt-8d", "Gt-CB1")

#For each sample...
message("Reading in functional annotation files")
for (strain in colnames(orthogroups)) {
  
  #Read in CSEP predictions
  CSEPs <- read.table(
    paste0("R:/GaeumannomycesGenomics/03_functional_annotation/", strain, "_CSEPs.tsv"),
    sep="\t", header=FALSE,
    col.names=c("Gene.ID", "Protein_ID", "PHI-base_entry",
                "Gene_name", "Pathogen_ID", "Organism", "Mutant_phenotype",
                "e-value", "bitscore", "pident"), fill=TRUE,
    comment.char=""
  )
  
  #Read in CAZyme predictions
  CAZymes <- read.table(
    paste0("S:/functional_annotation/001_run_dbcan/", strain, "/overview.txt"),
    sep="\t", header=TRUE,
    comment.char=""
  )
  
  #Read in AHRD tables
  AHRD <- read.table(
    paste0("S:/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/", strain, "/ahrd_output.csv"),
    sep="\t", header=TRUE, quote="",
    fill=TRUE, comment.char="#"
  )
  
  assign(paste0(strain, ".CSEPs"), CSEPs)
  assign(paste0(strain, ".CAZymes"), CAZymes)
  assign(paste0(strain, ".AHRD"), AHRD)
  
}

#Make dataframe with gene counts for each orthogroup
orthogroups.count <- orthogroups

for (i in 1:length(colnames(orthogroups.count))) {
  orthogroups.count[, i] <- sapply(strsplit(orthogroups[, i], " "), length)
}


## Assign AHRD results and biotype to orthogroups ##

#Combine all AHRD annotations
all.AHRD <- bind_rows(mget(ls(pattern=".AHRD"))) %>%
  mutate(Protein.Accession2=sub("-", "", Protein.Accession))

#Read in gene headers
gene.headers <- read.csv(
  "R:/GaeumannomycesGenomics/07_comparative_genomics/gene_list.txt",
  sep=" ", header=FALSE
  ) %>%
  rename(gene=V1, name=V2, note=V3, confidence=V4, representative=V5, biotype=V6) %>%
  mutate(name=sub("Name=", "", sub("-", "", name)),
         biotype=sub("biotype=", "", biotype))

#Make lists to capture annotations and biotype for HOGs
AHRD.list <- list()
biotype.list <- list()

#Create bar to show progress
progress.bar <- txtProgressBar(1, nrow(orthogroups), initial=0, char="=", style=3)

#For each HOG...
for (i in 1:nrow(orthogroups)) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, i)
  
  genes <- str_trim(unlist(str_split(orthogroups[i,], ",")))
  
  #Add AHRD annotation for all genes
  AHRD.list[[rownames(orthogroups[i,])]] <- 
    all.AHRD[which(all.AHRD$Protein.Accession2 %in% genes),]
  
  #Add biotype for all genes
  biotype.list[[rownames(orthogroups[i,])]] <- 
    gene.headers$biotype[which(gene.headers$name %in% genes)]
}
close(progress.bar)

#Label isoform variants with contradictory descriptions for appraisal
isoform.variants <- all.AHRD %>%
  mutate(variant=sub("\\..*$", "", Protein.Accession)) %>%
  group_by(variant) %>%
  mutate(n=n_distinct(Interpro.ID..Description.)) %>%
  ungroup() %>%
  arrange(desc(n))

#Proportion of variants with contradictory descriptions
isoform.variants %>%
  mutate(strain=sub("_EIv1.*$", "", Protein.Accession)) %>%
  group_by(strain, variant) %>%
  summarise(num.desc=max(n)) %>%
  ungroup() %>%
  group_by(strain) %>%
  summarise(genes=n_distinct(variant),
            contradictions=sum(num.desc > 1)) %>%
  mutate(proportion=contradictions/genes*100)


## Assign CSEP and CAZyme results to orthogroups ##

#Make dataframe for CSEP/CAZyme counts
count <- data.frame(
  matrix(0,
         ncol=ncol(orthogroups.count),
         nrow=nrow(orthogroups.count))
)
colnames(count) <- colnames(orthogroups.count)
rownames(count) <- rownames(orthogroups.count)

#Make list to capture genes not assigned to HOGs
missing.hogs <- list()

message("Counting number of CSEPs/CAZymes in each orthogroup:")

#For CSEPs and CAZymes in turn...
for (group in c("CSEP", "CAZyme")) {
  
  message(group, "s:")
  tmp.count <- count
  #Make list to capture metadata
  gene.list <- list()
  
  #For each sample...
  for (i in 1:length(colnames(tmp.count))) {
    
    #Print progress
    message((i - 1), "/", length(colnames(tmp.count)), " strains")
    
    #Retrieve the list of gene predictions
    genes <- get(paste0(colnames(tmp.count)[i], ".", group, "s"))
    
    #For each gene...
    for (j in 1:length(genes$Gene.ID)) {
      
      counter <- 1
      
      #Retrieve gene
      gene <- grep(sub("-", "", genes$Gene.ID[j]), orthogroups[, i])
      
      #Add to missing HOG list if not found
      if (length(gene) == 0) {
        
        missing.hogs[[group]][colnames(tmp.count)[i]][[counter]] <- sub("-", "", genes$Gene.ID[j])
        
        counter <- counter + 1
        
      } else {
        
        #Search for gene in corresponding column of orthogroups dataframe and add 1 to orthogroup count
        tmp.count[gene, i] <- tmp.count[gene, i] + 1
        
        #Get metadata
        gene.list[[rownames(tmp.count[gene,])]][[colnames(tmp.count)[i]]] <- genes[j,]
        
      }
      
    }
    
  }
  
  
  message(i, "/", length(colnames(tmp.count)), " strains")
  assign(paste0(group, ".count"), tmp.count)
  assign(paste0(group, ".list"), gene.list)
  
}


## Summarise orthogroup stats ##

#Read in sample metadata
metadata <- read.csv("R:/GaeumannomycesGenomics/05_phylogenomics/raxmlng/metadata.csv")

#Make dataframe to collect stats for orthogroups
orthogroups.stats <- data.frame(
  orthogroup=rownames(orthogroups.count),
  copy_number="multi",
  biotype=NA,
  CSEP=NA,
  CSEP_name=NA, 
  PHI.base_entry=NA,
  CAZyme=NA,
  EC=NA, 
  CAZyme_family=NA, 
  CAZyme_name=NA,
  PCWDE=NA,
  Gene.Ontology.Term=NA
)

#Add whether orthogroups are single copy
orthogroups.stats$copy_number[
  match(rownames(orthogroups.count[apply(orthogroups.count < 2, 1, all),]),
        orthogroups.stats$orthogroup)
] <- "single"


message("Assigning CSEP/CAZyme orthogroups as -only or -mixed:")

#Determine which orthogroups contain only CSEPs (SP-only) or both CSEPs and other genes (SP-mixed) 
CSEPs <- names(which(rowSums(CSEP.count) > 0))
CAZymes <- names(which(rowSums(CAZyme.count) > 0))

CSEP.type <- list()
CAZyme.type <- list()

#For each sample...
for (i in 1:length(colnames(orthogroups.count))) {
  
  #Print progress
  message((i - 1), "/", length(colnames(orthogroups.count)), " strains")
  
  for (CSEP in CSEPs) {
    if (orthogroups.count[CSEP,i] == CSEP.count[CSEP,i]) {
      CSEP.type[[CSEP]][i] <- "CSEP"
    } else {
      CSEP.type[[CSEP]][i] <- "mixed"
    }
  }
  
  for (CAZyme in CAZymes) {
    if (orthogroups.count[CAZyme,i] == CAZyme.count[CAZyme,i]) {
      CAZyme.type[[CAZyme]][i] <- "CAZyme"
    } else {
      CAZyme.type[[CAZyme]][i] <- "mixed"
    }
  }
  
}
message(i, "/", length(colnames(orthogroups.count)), " strains")

for (CSEP in CSEPs) {
  if(grepl("mixed", CSEP.type[CSEP]) == TRUE) {
    orthogroups.stats$CSEP[orthogroups.stats$orthogroup == CSEP] <- "mixed"
  } else {
    orthogroups.stats$CSEP[orthogroups.stats$orthogroup == CSEP] <- "CSEP"
  }
}

for (CAZyme in CAZymes) {
  if(grepl("mixed", CAZyme.type[CAZyme]) == TRUE) {
    orthogroups.stats$CAZyme[orthogroups.stats$orthogroup == CAZyme] <- "mixed"
  } else {
    orthogroups.stats$CAZyme[orthogroups.stats$orthogroup == CAZyme] <- "CAZyme"
  }
}


message("Adding AHRD GO terms:")

#Create bar to show progress
progress.bar <- txtProgressBar(1, nrow(orthogroups.stats), initial=0, char="=", style=3)

#For each orthogroup...
for (i in 1:nrow(orthogroups.stats)) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, i)
  
  go.terms <- unique(AHRD.list[[i]]$Gene.Ontology.Term)
  go.terms <- go.terms[go.terms != ""]
  
  if (length(go.terms) > 1) {
    
    orthogroups.stats$Gene.Ontology.Term[i] <-
      str_c(go.terms, collapse=" || ")
    
  } else if (length(go.terms) == 1) {
    
    orthogroups.stats$Gene.Ontology.Term[i] <-
      go.terms
    
  }
  
}
close(progress.bar)


message("Adding biotype:")

#Print orthogroups with multiple biotypes
#biotype.list[which(unlist(lapply(biotype.list, function(x) length(unique(x)))) > 1)]

#If any of the genes in an orthogroup are transposable element genes, conservatively assign orthogroup as such
orthogroups.stats$biotype[match(names(which(unlist(lapply(
  biotype.list, function(x) "transposable_element_gene" %in% x)))),
  orthogroups.stats$orthogroup)] <- "transposable_element_gene"


message("Adding CSEP/CAZyme annotations:")

#Create bar to show progress
progress.bar <- txtProgressBar(1, length(names(CAZyme.list)), initial=0, char="=", style=3)

#For each CAZyme...
for (i in 1:length(names(CAZyme.list))) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, i)
  
  #Get EC numbers
  EC.number <- unique(unlist(lapply(CAZyme.list[[i]], "[", "EC.")))
  
  #If gene is classified as more than one EC number...
  if (length(EC.number > 1)) {
    
    #Search for each EC number separately
    for (j in EC.number) {
      
      j <- sub("\\|.*", "", j)
      
      enzyme.list <- list()
      
      #Search ExplorEnz website for EC number
      explorenz <- read_html(paste0("https://www.enzyme-database.org/query.php?ec=", j))
      #Scrape accepted name
      enzyme.list[EC.number] <- explorenz %>%
        html_element("tr:nth-child(2) td+ td") %>% 
        html_text2()
      
    }
    
    #Summarise EC numbers/accepted names
    EC.number <- paste(EC.number, collapse=",")
    CAZyme.name <- paste(unique(unlist(enzyme.list)), collapse=",")
    
  } else {
    
    #Search ExplorEnz website for EC number
    explorenz <- read_html(paste0("https://www.enzyme-database.org/query.php?ec=", EC.number))
    #Scrape accepted name
    CAZyme.name <- explorenz %>%
      html_element("tr:nth-child(2) td+ td") %>% 
      html_text2()
    
  }
  
  #Add results to dataframe
  orthogroups.stats$EC[match(names(CAZyme.list)[i], orthogroups.stats$orthogroup)] <- EC.number
  orthogroups.stats$CAZyme_name[match(names(CAZyme.list)[i], orthogroups.stats$orthogroup)] <- CAZyme.name
  
  #Get CAZyme family
  CAZyme.family <- unique(unlist(lapply(CAZyme.list[[i]], "[", "DIAMOND")))
  
  #If gene is classified as more than one CAZyme family...
  if (length(CAZyme.family > 1)) {
    
    #Summarise families
    CAZyme.family <- paste(CAZyme.family, collapse=",")
    
  }
  
  #Add to dataframe
  orthogroups.stats$CAZyme_family[match(names(CAZyme.list)[i], orthogroups.stats$orthogroup)] <- CAZyme.family
  
}
close(progress.bar)

#Create bar to show progress
progress.bar <- txtProgressBar(1, length(names(CSEP.list)), initial=0, char="=", style=3)

#For each CSEP...
for (i in 1:length(names(CSEP.list))) {

  #Update progress bar
  setTxtProgressBar(progress.bar, i)

  #Get gene name
  CSEP.name <- unique(unlist(lapply(CSEP.list[[i]], "[", "Gene_name")))
  #Get PHI-base entry
  PHI.base <- unique(unlist(lapply(CSEP.list[[i]], "[", "PHI.base_entry")))

  #If there's an associated name...
  if (any(CSEP.name != "")) {

    #If gene is classified under more than one gene name...
    if (length(CSEP.name > 1)) {

      #Summarise names
      CSEP.name <- paste(CSEP.name[CSEP.name != ""], collapse=",")
      PHI.base <- paste(PHI.base[PHI.base != ""], collapse=",")

    }

    #Add results to dataframe
    orthogroups.stats$CSEP_name[match(names(CSEP.list)[i], orthogroups.stats$orthogroup)] <- CSEP.name
    orthogroups.stats$PHI.base_entry[match(names(CSEP.list)[i], orthogroups.stats$orthogroup)] <- PHI.base

  }

}
close(progress.bar)

message(paste0("Results saved in orthogroup-matrices-", Sys.Date(), ".RData"))
save(orthogroups,
     orthogroups.stats,
     orthogroups.count,
     CSEP.count,
     CAZyme.count,
     AHRD.list,
     file=paste0("orthogroup-matrices-", Sys.Date(), ".RData"))

write.csv(
  file=paste0(
    "R:/GaeumannomycesGenomics/07_comparative_genomics/orthogroup-stats-",
    Sys.Date(), ".csv"),
  orthogroups.stats, quote=FALSE)

#Write abundance matrix files for lifestyle PERMANOVA testing
write.csv(
  file="R:/GaeumannomycesGenomics/07_comparative_genomics/orthogroup-count.csv",
  orthogroups.count, quote=FALSE)
write.csv(
  file="R:/GaeumannomycesGenomics/07_comparative_genomics/CAZyme-count.csv",
  CAZyme.count, quote=FALSE)
write.csv(
  file="R:/GaeumannomycesGenomics/07_comparative_genomics/CSEP-count.csv",
  CSEP.count)
