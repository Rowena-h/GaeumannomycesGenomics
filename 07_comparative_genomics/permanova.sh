#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 1G                                # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

for group in orthogroup CSEP CAZyme BGC
do

	out_dir=lifestyle_permanova/${group}

	#Make output directory
	mkdir -p ${out_dir}

	singularity exec ~/programmes/gene_comparison/gene_comparison-230412.img python run_edited.py \
		-i ${group}_abundance_matrix.csv \
		-t species_tree_ingroup.tre \
		--colors nonpathogen:#009E73,pathogen:#696969 \
		-o ${out_dir}

done
