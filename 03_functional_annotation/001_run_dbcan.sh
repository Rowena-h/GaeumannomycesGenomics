#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 8                                    # number of cores
#SBATCH --mem 150000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

strain_file=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
protein_dir=../data/ncbi_data/proteins/gaeumannomyces
out_dir=../scratch/functional_annotation/001_run_dbcan/${strain_file}

#Make output directory
mkdir -p ${out_dir}

#seqkit v0.12.0
source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6

#Filter for representative isoforms
seqkit grep -r -n -p '.*representative=True.*' \
	../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/${strain_file}/*.release.gff3.pep.fasta \
	> ${protein_dir}/${strain_file}_EIv1.release.gff3.pep.repisoform.fasta

#HMMER v3.3
source package cd746169-51ea-44e5-b34c-4f9e80cbc470

#Run run_dbcan
singularity exec ~/programmes/run_dbcan/run_dbcan-3.0.1.img run_dbcan \
	${protein_dir}/${strain_file}_EIv1.release.gff3.pep.repisoform.fasta protein \
	--out_dir ${out_dir} --db_dir ~/programmes/run_dbcan/db
