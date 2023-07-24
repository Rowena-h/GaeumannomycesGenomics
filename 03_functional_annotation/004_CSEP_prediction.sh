#!/bin/bash
#SBATCH -p ei-medium				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 1					# number of cores
#SBATCH --mem 1G				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

strain_file=$(/usr/bin/awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/functional_annotation/004_CSEP_prediction

#Make output directory
mkdir -p $out_dir

proteins=../data/ncbi_data/proteins/gaeumannomyces/${strain_file}_EIv1.release.gff3.pep.repisoform.fasta
out_prefix=${strain_file}

for tool in deeploc deepsig effectorp1 effectorp2 effectorp3 phobius ps_scan signalp3 signalp4 signalp6 targetp tmhmm
do
	sbatch ${tool}.sh $proteins $out_prefix $out_dir
done

source blast-2.10

makeblastdb -in phi-base_current.fas -dbtype prot
