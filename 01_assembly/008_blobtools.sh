#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -n 4                                    # number of cores
#SBATCH --mem 100000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

source blobtools-1.0.1

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/008_blobtools/${strain_file}

#Make output directory
mkdir -p ${out_dir}

blobtools create	-i ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_ctg.fa \
                        -b ../scratch/007_readsmap/${strain_file}/${strain_file}_aln.sorted.bam \
                        -t ../scratch/006_blastn/${strain_file}/${strain_file}_hits \
                        --nodes ncbi_taxdump_090223/nodes.dmp \
                        --names ncbi_taxdump_090223/names.dmp \
                        -o ${out_dir}/${strain_file}

for rank in species family order
do

        blobtools view -i ${out_dir}/${strain_file}.blobDB.json -o ${out_dir}/${rank} -r ${rank}

        blobtools plot -i ${out_dir}/${strain_file}.blobDB.json -o ${out_dir}/${rank} -r ${rank}

done
