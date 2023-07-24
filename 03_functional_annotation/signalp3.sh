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

#seqkit v0.12.0
source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6

seqkit split $proteins -s 4000

source predector-1.2.7_CBG

rm ${out_dir}/${out_prefix}_signalp3.tsv
printf '# name\t!\tCmax\tpos\t?\tSprob\t?\n' > ${out_dir}/${out_prefix}_signalp3.tsv

for file in $(ls ${proteins}.split)
do
	signalp3 -type euk -method hmm -short ${proteins}.split/${file} | tail -n +3 >> ${out_dir}/${out_prefix}_signalp3.tsv
done
