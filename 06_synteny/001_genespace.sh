#!/bin/bash
#SBATCH -p ei-largemem                          # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 20                                   # number of cores
#SBATCH --mem 512G                              # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

#Format protein and gff files
for strain in Gt14LH10 Gt-19d1 Gt-23d Gt-4e Gt-8d Gh-1B17 Gh-2C17 Gt-3aA1 Gt-CB1
do
	mkdir -p ../results/hifiasm_assemblies/genespace/${strain}
	sed '/>/ s/ .*//; s/-//g' ../data/ncbi_data/proteins/gaeumannomyces/${strain}_EIv1.release.gff3.pep.repisoform.fasta > ../results/hifiasm_assemblies/genespace/${strain}/${strain}_pep_genespace.fasta
	sed '/^##/! s/-//g' ../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/${strain}/*_EIv1.release.gff3 > ../results/hifiasm_assemblies/genespace/${strain}/${strain}_genespace.gff3
done

#Run GENESPACE
singularity exec ~/programmes/genespace/230330_genespace.img Rscript genespace.R
