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

1. `sbatch -a 1-9 001_run_dbcan.sh` predicts CAZymes from protein sets using [run_dbcan](https://github.com/linnabrown/run_dbcan).
2. `sbatch -a 1-9 002_antismash.sh` predicts secondary metabolites from protein sets using [antiSMASH](https://github.com/antismash/antismash).
3. `sbatch -a 1-9 003_eggnog-mapper.sh` performs functional annotation of protein sets using [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper).
4. `sbatch -a 1-9 004_CSEP_prediction.sh` submits a suite of tools that feed into CSEP prediction: `deepsig.sh` ([DeepSig](https://github.com/BolognaBiocomp/deepsig)), `deeploc.sh` ([DeepLoc](https://services.healthtech.dtu.dk/services/DeepLoc-1.0/)), `effectorp1.sh` ([EffectorP v1](https://github.com/JanaSperschneider/EffectorP-1.0)), `effectorp2.sh` ([EffectorP v2](https://github.com/JanaSperschneider/EffectorP-2.0)), `effectorp3.sh` ([EffectorP v3](https://github.com/JanaSperschneider/EffectorP-3.0)), `phobius.sh` ([Phobius](https://phobius.sbc.su.se/)), [`ps_scan.sh`](https://github.com/ebi-pf-team/interproscan/blob/master/core/jms-implementation/support-mini-x86-32/bin/prosite/ps_scan.pl), `signalp3.sh` ([SignalP v3](https://services.healthtech.dtu.dk/services/SignalP-3.0/)), `signalp4.sh` ([SignalP v4.1](https://services.healthtech.dtu.dk/services/SignalP-4.1/)), `signalp6.sh` ([SignalP v6](https://services.healthtech.dtu.dk/services/SignalP-6.0/)), `targetp.sh` ([TargetP](https://services.healthtech.dtu.dk/services/TargetP-2.0/)), `tmhmm.sh` ([TMHMM](https://services.healthtech.dtu.dk/services/TMHMM-2.0/)).
6. `sbatch 005_CSEPfilter.sh` runs `CSEPfilter` to produce a list of CSEPs from the outputs of tools listed above.
7. `sbatch -a 1-9 006_CSEP_blastp.sh` searches CSEP sequences against the [PHI-base database](http://www.phi-base.org/) (requires `phi-base_current.fas` and `phi-base_current.csv` to be downloaded into this directory from [here](http://www.phi-base.org/downloadLink.htm)).

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
4. Scripts to plot figures: `plot_genespace.R`, `plot_read_coverage.R`
 
## 7 Comparative genomics

`cd 07_comparative_genomics` :file_folder:

1. `sbatch 001_orthogroup_assignment.sh` `orthogroup_assigner.R`
2. `sbatch 002_big-scape.sh`
3. `sbatch 003_lifestyle_test.sh` `lifestyle_v_phylogeny.R` `permanova.sh` `run_edited.py`
4. Scripts to plot figures: `plot_ideograms.R`, `plot_gene_content.R`

