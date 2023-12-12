#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 16                                   # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

strain_file=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
assembly=$(realpath ../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Analysis/${strain_file}/Reference/Genome/*_EI_v1.1_p_ctg.fa)
out_dir=../scratch/mitohifi/${strain_file}

#Make output directory
mkdir -p ${out_dir}

cd ${out_dir}

#Identify contigs corresponding to mitogenome
singularity exec /ei/software/cb/mitohifi/mitohifi-2.2/x86_64/mitohifi-2.2_CBG.sif mitohifi.py \
	-c ${assembly} \
	-f /ei/.project-scratch/d/d2c0bfb1-c37a-4211-9493-86b15d4e773e/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Reference/Mikado_Sordariomycetes/sordariomycetes_NC_072722.1.fasta \
	-g /ei/.project-scratch/d/d2c0bfb1-c37a-4211-9493-86b15d4e773e/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Reference/Mikado_Sordariomycetes/sordariomycetes_NC_072722.1.gb \
	-t ${SLURM_CPUS_PER_TASK} \
	-o 3 \
	-p 60 \
	--circular-size 250
