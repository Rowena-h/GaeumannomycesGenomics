#!/bin/bash
#SBATCH -p ei-long                              # queue
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

#Maps reads to assembly
minimap2 \
	-ax map-hifi -t ${SLURM_CPUS_PER_TASK} \
	../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_ctg.fa \
	../scratch/hifi-reads/${strain_file}.fastq.gz > ${out_dir}/${strain_file}_aln.sam

#Convert to sorted, indexed BAM
samtools view -S -b ${out_dir}/${sample}/${sample}_aln.sam > ${out_dir}/${sample}/${sample}_aln.bam
samtools sort ${out_dir}/${sample}/${sample}_aln.bam -o ${out_dir}/${sample}/${sample}_aln.sorted.bam
samtools index ${out_dir}/${sample}/${sample}_aln.sorted.bam
