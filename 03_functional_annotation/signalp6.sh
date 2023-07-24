#!/bin/bash
#SBATCH -p ei-medium				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 1					# number of cores
#SBATCH --mem 10G				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

source predector-1.2.7_CBG

proteins=$1
out_prefix=$2
out_dir=$3

signalp6 \
    --fastafile $proteins \
    --output_dir ${out_dir}/${out_prefix}_signalp6 \
    --format none \
    --organism eukarya \
    --mode fast \
    --bsize 64

tail +2 ${out_dir}/${out_prefix}_signalp6/prediction_results.txt > ${out_dir}/${out_prefix}_signalp6.tsv
