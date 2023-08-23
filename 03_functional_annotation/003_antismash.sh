#!/bin/bash
#SBATCH -p ei-medium                     	# queue
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 6                            	# number of cores
#SBATCH --mem 150000	                     	# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail            	# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk 	# send-to address

#R v3.6
source package e03556a3-2fc2-4dc0-842d-fda34158940e

strain_file=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)

annotation_dir=../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/${strain_file}
out_dir=../scratch/functional_annotation/002_antismash/${strain_file}

#Make output directory
mkdir -p ${out_dir}

#Remove non-representative variants
Rscript remove_variants.r ${annotation_dir}/${strain_file}_EIv1.release.gff3 ${strain_file}_edited ../scratch/functional_annotation/002_antismash/

singularity exec ~/programmes/antismash/antismash-6.1.1.img antismash \
	../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Analysis/${strain_file}/Reference/Genome/*ctg.fa \
	--output-dir ${out_dir} \
	--taxon fungi \
	--cassis \
	--cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go \
	--output-basename ${strain_file}_antismash \
	--genefinding-gff3 ../scratch/functional_annotation/002_antismash/${strain_file}_edited.gff3 \
	--allow-long-headers \
	-v \
	-c ${SLURM_CPUS_PER_TASK}
