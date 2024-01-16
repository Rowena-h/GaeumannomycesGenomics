#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 1GB	                        # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

strain=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
new_strain=$(awk '{print $4}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
locus_tag=$(awk '{print $5}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
biosample=$(awk '{print $6}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
assembly=$(ls ../data_upload/ncbi_assemblies/*${new_strain}*ctg.fa)
out_dir=../data_upload/ncbi_annotation

#Make output directory
mkdir -p ${out_dir}

#Make NCBI submission template for strain
sed "s/STRAIN/${biosample}/" ncbi_template.sbt > ${out_dir}/${new_strain}_template.sbt

~/programmes/table2asn/linux64.table2asn \
        -t ${out_dir}/${new_strain}_template.sbt \
        -locus-tag-prefix ${locus_tag} \
        -euk \
        -gaps-unknown 100 \
        -i ${assembly} \
        -f ../results/ncbi_submission/${strain}/${strain}_merged_variants.gff3 \
        -o ${out_dir}/${new_strain}.sqn \
        -M n -V b -Z
