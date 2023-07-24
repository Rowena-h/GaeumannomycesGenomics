#!/bin/bash
#SBATCH -p ei-short				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 1					# number of cores
#SBATCH --mem 1G				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

for strain in Gt-CB1 Gt-3aA1 Gh-2C17 Gh-1B17 Gt-23d Gt-8d Gt-19d1 Gt-4e Gt14LH10
do
	./CSEPfilter ${strain} ../scratch/functional_annotation/004_CSEP_prediction
done
