#!/bin/bash
#SBATCH -p ei-medium                     	# queue
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 6                            	# number of cores
#SBATCH --mem 150000	                     	# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail            	# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk 	# send-to address

out_dir=../scratch/functional_annotation/002.5_big-scape/gaeumannomyces

#Make output directory
mkdir -p ${out_dir}

strains=$(awk '{printf "%s%s",sep,$1; sep=" "} END{print ""}' ../strains)

singularity exec ~/programmes/big-scape/big-scape.img python /usr/src/BiG-SCAPE/bigscape.py \
	-i ../scratch/functional_annotation/002_antismash \
	--include_gbk_str [ ${strains} ] \
	-o ${out_dir} \
	-c ${SLURM_CPUS_PER_TASK}
