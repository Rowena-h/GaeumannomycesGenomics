#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 16                                   # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

source lmod-6.1
ml kat/2.3.4

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/003_kat_comp/${strain_file}

#Make output directory
mkdir -p ${out_dir}

kat comp \
	-m 63 -t ${SLURM_CPUS_PER_TASK} \
	-o ${out_dir}/${strain_file}_hifiasm_vs_hifi \
	-H 1000000000 <(gunzip -c ../scratch/hifi-reads/*${strain_file}_hifi.fastq.gz) ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_ctg.fa
