#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 16                                   # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

source lmod-6.1
ml kat/2.3.4

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
orig_strain=$(echo ${strain_file} | sed 's/_.*$//')
out_dir=../scratch/011_kat_comp/${strain_file}
asm_file=../scratch/010_filter_lowcov/${strain_file}/${strain_file}.asm.bp.p_ctg.hifi_sect_filtered.fa
reads_file=../scratch/hifi-reads/*${strain_file}_hifi.fastq.gz

#Make output directory
mkdir -p ${out_dir}

kat comp \
	-m 63 -t ${SLURM_CPUS_PER_TASK} \
	-o ${out_dir}/${strain_file}_filteredhifiasm_vs_hifi \
	-H 1000000000 <(gunzip -c ${reads_file}) ${asm_file}

kat plot spectra-cn ${out_dir}/${strain_file}_filteredhifiasm_vs_hifi-main.mx \
	-x 300 \
	-o ${out_dir}/${strain_file}_filteredhifiasm_vs_hifi-main.mx.spectra-cn.pdf
