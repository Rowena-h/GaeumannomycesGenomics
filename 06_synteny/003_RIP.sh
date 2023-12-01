#!/bin/bash
#SBATCH -p ei-medium                     	# queue
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 1                            	# number of cores
#SBATCH --mem 150000	                     	# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail            	# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk 	# send-to address

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
strain=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/RIP/${strain_file}

#Make output directory
mkdir -p ${out_dir}

#seqkit
source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6

awk -v var="${strain}" '{if($1 == var) {print $2}}' ../06_synteny/pseudochromosomes.tsv | sed 's/.*_//' > ${strain_file}_tmp

seqkit grep -nrf ${strain_file}_tmp ../scratch/010_filter_lowcov/${strain_file}/${strain_file}.asm.bp.p_ctg.hifi_sect_filtered.fa > ${out_dir}/${strain_file}_chr.fa

rm ${strain_file}_tmp

singularity exec ~/programmes/BioPerl/RIP.img perl /bin/RIP_index_calculation.pl \
	-i ${out_dir}/${strain_file}_chr.fa \
	-r bed \
	-t CRI \
	-w 500 \
	-s 500 > ${out_dir}/${strain_file}_RIP.bed

~/scripts/RIP_genome_summary.pl ${out_dir}/${strain_file}_RIP.bed > ${out_dir}/${strain_file}_RIP_stats.txt
