#!/bin/bash
#SBATCH -p ei-medium                            # Use normal partition (queue) for now.
#SBATCH -N 1                                    # number of nodes
#SBATCH -n 1                                    # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/014_tidk/${strain_file}

#Make output directory
mkdir -p ${out_dir}

asm_file=../scratch/010_filter_lowcov/${strain_file}/${strain_file}.asm.bp.p_ctg.hifi_sect_filtered.fa
reads_file=../scratch/hifi-reads/*${strain_file}_hifi.fastq.gz

singularity exec ~/programmes/tidk/tidk.img tidk search \
	-s TTAGGG \
	-o ${strain_file} \
	-d ${out_dir} \
	-e tsv \
	${asm_file}

singularity exec ~/programmes/tidk/tidk.img tidk plot \
	-t ${out_dir}/${strain_file}_telomeric_repeat_windows.tsv \
	-o ${out_dir}/${strain_file}

