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
out_dir=../scratch/009_kat_sect/${strain_file}

#Use contaminant filtered unitigs if they exist, otherwise the contigs
if [[ -f ../scratch/009_remove_contam/${strain_file}/${strain_file}.asm.bp.p_utg.fa ]]
then

	asm_file=../scratch/009_remove_contam/${strain_file}/${strain_file}.asm.bp.p_utg.fa
	reads_file=$(awk '{print $3}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)

else

	asm_file=../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_ctg.fa
	reads_file=../scratch/hifi-reads/*${strain_file}_hifi.fastq.gz

fi

#Make output directory
mkdir -p ${out_dir}

kat sect \
	-m 63 -t ${SLURM_CPUS_PER_TASK} -n -o ${out_dir}/${strain_file}_hifi_vs_assm_sect \
	${asm_file} <(gunzip -c ${reads_file})
