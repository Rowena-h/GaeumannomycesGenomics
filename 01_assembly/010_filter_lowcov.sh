#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
orig_strain=$(echo ${strain_file} | sed 's/_.*$//')
out_dir=../scratch/010_filter_lowcov/${strain_file}

#Make output directory
mkdir -p ${out_dir}

#Use contaminant filtered unitigs if they exist, otherwise the contigs
if [[ -f ../scratch/009_remove_contam/${strain_file}/${strain_file}.asm.bp.p_utg.fa ]]
then
	
	asm_file=../scratch/009_remove_contam/${strain_file}/${strain_file}.asm.bp.p_utg.fa

	#Filter low coverage contigs according to KAT
	singularity exec ~/singularity_imgs/Bioconductor.img Rscript scripts/low_cov_deleter2.R \
        	${asm_file} \
        	../scratch/002_kat_hist/${orig_strain}/${orig_strain}_63mers \
        	../scratch/009_kat_sect/${strain_file}/${strain_file}_hifi_vs_assm_sect-stats.tsv \
        	${out_dir}/${strain_file}.asm.bp.p_utg.hifi_sect_filtered.fa

else

        asm_file=../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_ctg.fa

	#Filter low coverage contigs according to KAT
        singularity exec ~/singularity_imgs/Bioconductor.img Rscript scripts/low_cov_deleter2.R \
                ${asm_file} \
                ../scratch/002_kat_hist/${orig_strain}/${orig_strain}_63mers \
                ../scratch/009_kat_sect/${strain_file}/${strain_file}_hifi_vs_assm_sect-stats.tsv \
                ${out_dir}/${strain_file}.asm.bp.p_ctg.hifi_sect_filtered.fa

fi
