#!/bin/bash
#SBATCH -p ei-medium				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 1					# number of cores
#SBATCH --mem 1G				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

proteins=$1
out_prefix=$2
out_dir=$3

~/programmes/ps_scan/ps_scan.pl -p PS00014 -o scan -d ~/programmes/ps_scan/prosite.dat $proteins > ${out_dir}/${out_prefix}_ps_scan.tsv
