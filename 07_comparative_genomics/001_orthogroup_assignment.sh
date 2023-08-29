#!/bin/bash
#SBATCH -p ei-medium                     	# queue
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 1                            	# number of cores
#SBATCH --mem 150000	                     	# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail            	# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk 	# send-to address

grep -h ">" ../data/ncbi_data/proteins/gaeumannomyces/*.fasta | sed 's/^>//' > gene_list.txt

singularity exec ~/programmes/gene_comparison/gene_comparison-230412.img Rscript orthogroup_assigner.r
