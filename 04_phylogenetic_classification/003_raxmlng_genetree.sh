#!/bin/bash
#SBATCH -p ei-long                              # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 2                                    # number of cores
#SBATCH --mem 1GB                               # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

for marker in $(cat markers)
do

	if [[ ${marker} = avenacinase ]]
	then
		model=JTT+I+G4
	else
		model=GTR+G
	fi

	aln=alignments/${marker}_genetree_aln.fasta

	if [[ -f "alignments/${marker}_genetree.fasta" ]]
	then
		#Convert file for RAxML-NG
		singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
			--parse \
	        	--msa ${aln} \
			--model ${model} \
	        	--prefix raxmlng/gaeumannomyces_${marker}

		if [[ ${marker} = ITS2 ]]
		then
			msa=raxmlng/gaeumannomyces_${marker}.raxml.reduced.phy
		else
			msa=${aln}
		fi		

		#Run ML tree search and bootstrapping for <=1000 iterations
		singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
			--all \
			--msa ${msa} \
			--model ${model} \
			--prefix raxmlng/gaeumannomyces_${marker} \
			--seed 2 \
			--threads auto{${SLURM_CPUS_PER_TASK}} \
			--bs-trees autoMRE{1000}			

		#Check convergence
		singularity exec ~/programmes/raxml-ng/raxml-ng-1.1.0.img raxml-ng \
			--bsconverge \
		        --bs-trees raxmlng/gaeumannomyces_${marker}.raxml.bootstraps \
		        --prefix raxmlng/gaeumannomyces_${marker}_convergence_test \
		        --seed 2
	fi

done
