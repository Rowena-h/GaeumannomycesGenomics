#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 1GB	                        # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

strain=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
organism=$(awk '{print $3}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
new_strain=$(awk '{print $4}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
dir=$(ls ../results/hifiasm_assemblies/*/v1.1/sequences/*${strain}* | sed -n 1p | sed 's|\(.*\)/.*|\1|')

#Merge into chromosomes
singularity exec ~/programmes/agptools/agptools.img agptools assemble \
	${dir}/*${strain}_EI_v1.1_p_ctg.fa \
	${dir}/${strain}_EI_v1.1_p_ctg_curated.agp > ${dir}/${strain}_merged.fa

#Rename headers
for chr in $(grep ">chr" ${dir}/${strain}_merged.fa | sed 's/>//')
do
	num=$(echo $chr | sed 's/chr//')
	
	sed -i "s/${chr}/${chr} [organism=${organism}] [strain=${new_strain}] [location=chromosome] [chromosome=${num}]/" ${dir}/${strain}_merged.fa
done	

sed -i "/>ptg/ s/$/ [organism=${organism}] [strain=${new_strain}]/" ${dir}/${strain}_merged.fa
