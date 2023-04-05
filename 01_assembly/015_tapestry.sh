#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 5                                    # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/015_tapestry/${strain_file}

#Make output directory
mkdir -p ${out_dir}

#Use contaminant, coverage filtered unitigs if they exist, otherwise the contigs
if [[ -f ../scratch/010_filter_lowcov/${strain_file}/${strain_file}.asm.bp.p_utg.hifi_sect_filtered.fa ]]
then

        asm_file=../scratch/010_filter_lowcov/${strain_file}/${strain_file}.asm.bp.p_utg.hifi_sect_filtered.fa
        reads_file=$(awk '{print $3}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)

else

        asm_file=../scratch/010_filter_lowcov/${strain_file}/${strain_file}.asm.bp.p_ctg.hifi_sect_filtered.fa
        reads_file=../scratch/hifi-reads/*${strain_file}_hifi.fastq.gz

fi

singularity exec ~/programmes/tapestry/tapestry.img weave \
	-a ${asm_file} \
	-r ../scratch/hifi-reads/*${strain_file}_hifi.fastq.gz \
	-t TTAGGG \
	-o ${out_dir} \
	-c ${SLURM_CPUS_PER_TASK}
