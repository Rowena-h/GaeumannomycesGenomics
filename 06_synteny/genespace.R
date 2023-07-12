library(GENESPACE)

wd.gaeumannomyces <- "/ei/.project-scratch/d/d2c0bfb1-c37a-4211-9493-86b15d4e773e/genespace_gaeumannomyces_230426"
wd.gt <- "/ei/.project-scratch/d/d2c0bfb1-c37a-4211-9493-86b15d4e773e/genespace_gt_230426"

genomeRepo <- "/ei/projects/d/d2c0bfb1-c37a-4211-9493-86b15d4e773e/GaeumannomycesGenomics/results/hifiasm_assemblies/genespace"

strains.gt <- 
  c("Magpoae",
    "Gt19d1",
    "Gt8d",
    "Gt23d",
    "Gt4e",
    "Gt14LH10")

files.gt <- 
  c("Magpoae",
    "Gt-19d1",
    "Gt-8d",
    "Gt-23d",
    "Gt-4e",
    "Gt14LH10")

strains.gaeumannomyces <- 
  c("Magpoae",
    "Gh1B17",
    "Gh2C17",
    "GtCB1",
    "Gt3aA1",
    "Gt19d1",
    "Gt8d",
    "Gt23d",
    "Gt4e",
    "Gt14LH10")

files.gaeumannomyces <- 
  c("Magpoae",
    "Gh-1B17",
    "Gh-2C17",
    "Gt-CB1",
    "Gt-3aA1",
    "Gt-19d1",
    "Gt-8d",
    "Gt-23d",
    "Gt-4e",
    "Gt14LH10")

for (i in c("gaeumannomyces", "gt")) {
  
  wd <- get(paste0("wd.", i))
  strains <- get(paste0("strains.", i))
  files <- get(paste0("files.", i))
  
  parsedPaths <- parse_annotations(
    rawGenomeRepo=genomeRepo,
    genomeDirs=files,
    genomeIDs=strains,
    gffString="_genespace.gff3",
    faString="_pep_genespace.fasta",
    headerEntryIndex=1,
    gffIdColumn="ID",
    genespaceWd=wd)
  
  gpar <- init_genespace(
    wd=wd,
    genomeIDs=strains,
    outgroup="Magpoae",
    path2mcscanx="/opt/MCScanX/",
    path2diamond="/opt/diamond"
  )
  
  out <- run_genespace(gpar)
  
}
