#!/bin/bash
#SBATCH -p ei-short                      	# queue
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 1                            	# number of cores
#SBATCH --mem 1000	                     	# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail            	# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk 	# send-to address

out_dir=alignments

#Format trimmed files
#rename N0. N0_ ${out_dir}/orthogroups/*_trimmed.fa

#Concatenate
python ~/scripts/AMAS.py concat -f fasta \
				-d aa \
				-i ${out_dir}/orthogroups/*_aln_trimmed.fa \
				-p ${out_dir}/gaeumannomyces_partition.txt \
				-t ${out_dir}/gaeumannomyces_concat.fa \
				-u fasta

#Add gene models
sed -i 's/^/JTT+I+G4, /' ${out_dir}/gaeumannomyces_partition.txt
