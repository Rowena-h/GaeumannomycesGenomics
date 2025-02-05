#!/bin/bash
#SBATCH -p ei-short                     	# queue
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 1                            	# number of cores
#SBATCH --mem 150000	                     	# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail            	# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk 	# send-to address

for bgc in $(ls bgc_gbks/)
do
	singularity exec ~/programmes/clinker/clinker.img clinker bgc_gbks/${bgc}/*.gbk \
		-p bgc_gbks/${bgc}/${bgc}.html
done
