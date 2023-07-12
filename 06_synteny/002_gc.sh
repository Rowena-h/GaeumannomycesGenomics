#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 1G                                # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

#samtools v1.15.1
source package 638df626-d658-40aa-80e5-14a275b7464b
#bedtools v2.29.2
source package 6394519c-541f-479c-b064-dd0b912eac04

out_dir=../scratch/gc

#Make output directory
mkdir -p ${out_dir}

#Calculate GC content in 1000 bp windows
for strain in Gt14LH10 Gt-19d1 Gt-23d Gt-4e Gt-8d Gh-1B17 Gh-2C17 Gt-3aA1 Gt-CB1
do
	samtools faidx ../results/hifiasm_assemblies/*${strain}/v1.1/sequences/*${strain}_EI_v1.1_p_ctg.fa
	bedtools makewindows 	-g ../results/hifiasm_assemblies/*${strain}/v1.1/sequences/*${strain}_EI_v1.1_p_ctg.fa.fai \
				-w 1000 > ${out_dir}/${strain}_windows_1000.bed
	bedtools nuc 	-fi ../results/hifiasm_assemblies/*${strain}/v1.1/sequences/*${strain}_EI_v1.1_p_ctg.fa  \
			-bed ${out_dir}/${strain}_windows_1000.bed > ${out_dir}/${strain}_gc.tsv
done
