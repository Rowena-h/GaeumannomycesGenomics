#!/bin/bash
#SBATCH -p ei-long                              # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 16                                   # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

#BUSCO v2.5.1
source package 97a7f391-0f5c-4bdb-b951-14d0d5ee4576

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/005_busco_asco/${strain_file} 

#Make output directory
mkdir -p ${out_dir}

#Run BUSCO
busco 	-m genome \
	-i ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_ctg.fa \
	-o hifi_assm_${strain_file}_busco_asco \
	-l ascomycota --out_path ${out_dir} \
	-c ${SLURM_CPUS_PER_TASK} -f \
	--offline --download_path ../data/busco_database/v5/data
