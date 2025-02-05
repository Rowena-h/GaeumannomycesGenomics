#!/bin/bash
#SBATCH -p ei-short                      	# queue
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 1                            	# number of cores
#SBATCH --mem 1000	                     	# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail            	# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk 	# send-to address

out_dir=tyr_genetree
marker=tyr

source mafft-7.271

#Create alignment
mafft ${out_dir}/${marker}_genetree.fasta > ${out_dir}/${marker}_genetree_aln.fasta

#Trim alignment
singularity exec ~/programmes/trimAl/trimAl.img trimal \
	-in ${out_dir}/${marker}_genetree_aln.fasta \
	-fasta -gappyout > ${out_dir}/${marker}_genetree_aln_trim.fasta

#Convert leading or trailing gaps to ?s
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${out_dir}/${marker}_genetree_aln_trim.fasta | \
	awk -F '[^-](.*[^-]|$)' '{s=$0; h=gsub(/./,"?",$1); t=gsub(/./,"?",$2); print $1 substr(s,h+1, length(s)-h-t) $2}' > tmp.fa && mv tmp.fa ${out_dir}/${marker}_genetree_aln_trim.fasta
