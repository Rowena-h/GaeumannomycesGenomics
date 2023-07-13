# GaeumannomycesGenomics

README IN PROGRESS...
Insert schematic

> Paper citation

The pipeline was written for and run on Norwich BioScience Institutes' HPC cluster which uses the SLURM batch-queue system. This means that many of the bash scripts (`.sh` file endings) specify core allocation, run times and memory usage allocation that may need to be adapted for different platforms.

---

## 1 *De novo* genome assembly

`cd 01_assembly` :file_folder:

### Checking reads

1. `sbatch -a 1-9 001_kat_hist.sh` produces kmer histograms for the raw HiFi reads using [KAT](https://github.com/TGAC/KAT).

### Assembly

2. `sbatch -a 1-9 002_hifiasm.sh` assembles HiFi reads using [hifiasm](https://github.com/chhylp123/hifiasm).
3. `sbatch -a 1-9 003_kat_comp.sh` checks for content correctness of assemblies with respect to the input HiFi reads using KAT.

### Assessment

4. `sbatch 004_quast.sh` produces assembly contiguity statistics using [QUAST](https://github.com/ablab/quast).
5. `sbatch -a 1-9 005_busco_asco.sh` produces gene set completeness statistics using [BUSCO](https://busco.ezlab.org/).
6. `sbatch -a 1-9 006_blastn.sh` searches the assembly against the NCBI nucleotide database using [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) to produce input for BlobTools.
7. `sbatch -a 1-9 007_readsmap.sh` maps HiFi reads to the assembly using [minimap2](https://github.com/lh3/minimap2) to produce input for BlobTools.
8. `sbatch -a 1-9 008_blobtools.sh` runs a contamination check using [BlobTools](https://github.com/DRL/blobtools).

### Filtering

9. `sbatch -a 1-9 009_kat_sect.sh` estimates assembly coverage levels using KAT.
10. `sbatch -a 1-9 010_filter_lowcov.sh` filters out small, low coverage sequences using a custom script, `scripts/low_cov_deleter.R`.

### Final assessment

11. `sbatch -a 1-9 011_kat_comp.sh` reruns KAT comp to check final content correctness of the filtered assemblies.
13. `sbatch -a 1-9 013_busco_asco.sh` reruns BUSCO on the filtered assemblies.
14. `sbatch 014_quast.sh` reruns QUAST on the filtered assemblies.

## 2 Structural annotation

## 3 Functional annotation

`cd 03_functional_annotation` :file_folder:

1. `sbatch -a 1-9 001_run_dbcan.sh` estimates CAZymes from predicted protein sets using [run_dbcan](https://github.com/linnabrown/run_dbcan).

## 4 Phylogenetic classification

`cd 04_phylogenetic_classification` :file_folder:

This folder contains a file - `markers` - listing the genetic markers selected for building the trees.

### Gene extraction

1. `sbatch 001_genepull.sh` uses [GenePull](https://github.com/Rowena-h/MiscGenomicsTools/tree/main/GenePull) to extract selected genetic markers from filtered assemblies of each strain. Requires fasta files containing a single example sequence from a closely related taxon for each genetic marker being extracted.

### Alignment and ML tree building

2. `sbatch 002_align_multicopy.sh` makes gene alignments for gdo and ITS2 for distinguishing *Gaeumannomyces* genetic groups using [MAFFT](https://github.com/GSLBiotech/mafft).
3. `sbatch 003_raxmlng_genetree.sh` builds ML gene trees for each marker using [RAxML-NG](https://github.com/amkozlov/raxml-ng) with bootstrapping until convergence or up to 1,000 replicates (whichever first).
4. `sbatch 004_align_singlecopy.sh` uses MAFFT to align single-copy markers and a single copy of the multi-copy markers from the previous alignments. Gene alignments are manually checked with [AliView](https://github.com/AliView/AliView), and `file_prep.sh` contains example one-liners for formatting sequence headers in each of the gene alignment fasta files so that they are identical across different genes (i.e. removing GenBank accessions; removing misc text after taxon names/vouchers; replacing spaces with underscores etc).
5. `sbatch 005_concat.sh` concatenates gene alignments using [AMAS](https://github.com/marekborowiec/AMAS).
6. `sbatch 006_raxmlng_speciestree.sh` builds ML species trees using RAXML-NG with bootstrapping until convergence or up to 1,000 replicates (whichever first).

## 5 Phylogenomics

`cd 05_phylogenomics` :file_folder:

1. `sbatch 001_orthofinder.sh` infers phylogenetic hierarchical orthogroups (HOGs) using [OrthoFinder](https://github.com/davidemms/OrthoFinder).
2. `sbatch 002_align_singlecopy.sh` submits batch array jobs to align all single-copy HOGs using MAFFT and trim using [trimAl](http://trimal.cgenomics.org/).
3. `sbatch 003_concat.sh` concatenates single-copy HOG alignments using AMAS.
4. `sbatch 004_raxmlng.sh` builds genome-scale ML species trees using RAXML-NG with bootstrapping until convergence or up to 1,000 replicates (whichever first).
5. Script to plot figure: `plot_trees.R`

## 6 Synteny and structure

`cd 06_synteny` :file_folder:

1. `sbatch 001_genespace.sh` formats protein and gff3 files and submits `genespace.R` to infer synteny between strains using [GENESPACE](https://github.com/jtlovell/GENESPACE).
2. `sbatch 002_gc.sh` calculates GC content in 1,000 bp windows across each genome using [bedtools](https://github.com/arq5x/bedtools2).
3. `sbatch 003_contigs2pseudochromosomes.sh` replaces fragment names according to pseudochromosomes inferred from GENESPACE, as recorded in `pseudochromosomes.tsv`.
4. Scripts to plot figures: `plot_ideograms.R`, `plot_genespace.R`, `plot_read_coverage.R`
