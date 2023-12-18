#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 5                                    # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

source minimap2-2.21
source samtools-1.13

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/007_readsmap/${strain_file}

#Make output directory
mkdir -p ${out_dir}

#Remove secondary alignments
samtools view -h ${out_dir}/${strain_file}_aln.sorted.bam | awk '$10 != "*"' | samtools sort > ${out_dir}/${strain_file}_aln.sorted.filtered.bam
samtools index ${out_dir}/${strain_file}_aln.sorted.filtered.bam
