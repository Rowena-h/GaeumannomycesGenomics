#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 8                                    # number of cores
#SBATCH --mem 150000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

source eggnog-mapper-2.1.9_CBG

strain_file=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/functional_annotation/003_eggnog-mapper/${strain_file}
proteins=$(ls ../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/${strain_file}/*.release.gff3.pep.fasta)

#Make output directory
mkdir -p ${out_dir}

emapper.py \
	--cpu ${SLURM_CPUS_PER_TASK} \
	--itype proteins \
	-m diamond \
	--pfam_realign none \
	--go_evidence non-electronic \
	--report_orthologs \
	-d /ei/cb/common/Databases/eggnog/eggnog_proteins.dmnd \
	--report_no_hits --override \
	-i ${proteins} \
	-o ${out_dir}/${strain_file}
