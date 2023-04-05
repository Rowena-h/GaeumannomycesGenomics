#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 2                                    # number of cores
#SBATCH --mem 10000                             # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

source raxml-ng-0.8.0

for marker in $(cat markers)
do

	if [[ -f "alignments/${marker}_genetree.fasta" ]]
	then
		#Convert file for RAxML-NG
		raxml-ng --parse \
	        	 --msa alignments/${marker}_genetree_aln.fasta \
			 --model GTR+G \
	        	 --prefix raxmlng/gaeumannomyces_${marker}

		if [[ ${marker} = ITS2 ]]
		then		

			#Run ML tree search and bootstrapping for <=1000 iterations
			raxml-ng --all \
				 --msa raxmlng/gaeumannomyces_${marker}.raxml.reduced.phy \
				 --model GTR+G \
			         --prefix raxmlng/gaeumannomyces_${marker} \
			         --seed 2 \
			         --threads ${SLURM_CPUS_PER_TASK} \
		        	 --bs-trees autoMRE{1000}

		else

			#Run ML tree search and bootstrapping for <=1000 iterations
	                raxml-ng --all \
                        	 --msa alignments/${marker}_genetree_aln.fasta \
                        	 --model GTR+G \
                        	 --prefix raxmlng/gaeumannomyces_${marker} \
                        	 --seed 2 \
                        	 --threads ${SLURM_CPUS_PER_TASK} \
                        	 --bs-trees autoMRE{1000}			

		fi

		#Check convergence
		raxml-ng --bsconverge \
		         --bs-trees raxmlng/gaeumannomyces_${marker}.raxml.bootstraps \
		         --prefix raxmlng/gaeumannomyces_${marker}_convergence_test \
		         --seed 2
	fi

done
