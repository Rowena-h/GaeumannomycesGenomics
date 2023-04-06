#!/bin/bash
#SBATCH -p ei-long				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 8					# number of cores
#SBATCH --mem 150000				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

source gcc-6.2.0
source zlib-1.2.11

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/hifiasm-assemblies/${strain_file}

#Make output directory
mkdir -p ${out_dir}

#Assemble reads
~/programmes/hifiasm/hifiasm \
	-o ${out_dir}/${strain_file}.asm \
	-t ${SLURM_CPUS_PER_TASK} \
	-l0 ../scratch/hifi-reads/${strain_file}.fastq.gz

#Convert hifiasm GFA files to fasta files
awk '/^S/{print ">"$2;print $3}' ${out_dir}/${strain_file}.asm.bp.p_ctg.gfa > ${out_dir}/${strain_file}.asm.bp.p_ctg.fa
