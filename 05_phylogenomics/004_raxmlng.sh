#!/bin/bash
#SBATCH -p ei-long                      	# Use normal partition (queue) for now.
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 48                            	# number of cores per task
#SBATCH --mem 100000	                     	# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail            	# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk 	# send-to address

alignment_dir=alignments
out_dir=raxmlng

#Make output directory
mkdir -p ${out_dir}

#Convert file for RAxML-NG
singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
	--parse \
        --msa ${alignment_dir}/gaeumannomyces_concat.fa \
        --model ${alignment_dir}/gaeumannomyces_partition.txt \
        --prefix ${out_dir}/gaeumannomyces_concat

#Run ML tree search and bootstrapping for <=1000 iterations
singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
	--all \
	--msa ${alignment_dir}/gaeumannomyces_concat.fa \
        --model ${alignment_dir}/gaeumannomyces_partition.txt \
        --prefix ${out_dir}/gaeumannomyces_concat \
        --seed 2 \
        --threads ${SLURM_CPUS_PER_TASK} \
        --bs-trees autoMRE{1000}

#Check convergence
singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
	--bsconverge \
        --bs-trees ${out_dir}/gaeumannomyces_concat.raxml.bootstraps \
        --prefix ${out_dir}/gaeumannomyces_concat_convergence_test \
        --seed 2
