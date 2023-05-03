#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 10                                   # number of cores
#SBATCH --mem 400000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

source blast-2.10

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
out_dir=../scratch/006_blastn/${strain_file}

#Make output directory
mkdir -p ${out_dir}

#BLAST assembly against nt database
blastn \
    -task megablast \
    -db /ei/public/databases/blast/ncbi/nt_20210521/nt \
    -query ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_ctg.fa \
    -out ${out_dir}/${strain_file}_hits \
    -outfmt "6 qseqid staxids bitscore std" \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -evalue 1e-25
