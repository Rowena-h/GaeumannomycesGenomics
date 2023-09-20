#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 4                                    # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

#QUAST v5.0.2
source package 65df873c-d601-44ba-ac61-64644b55dfbb

asm_files=$(awk '{print $2}' ../strains | sed 's/$/\.asm\.bp\.p_ctg\.hifi_sect_filtered\.fa/' | sed 's|^|\.\./scratch/010_filter_lowcov/\*/|')
strains=$(awk '{print $1}' ../strains)

#Run QUAST for contiguity statistics
quast	${asm_files} \
	-l $(echo ${strains}| sed 's/ /,/g') \
	-o ../scratch/014_quast/gaeumannomyces_final_quast_results \
	-t ${SLURM_CPUS_PER_TASK} --fungus
