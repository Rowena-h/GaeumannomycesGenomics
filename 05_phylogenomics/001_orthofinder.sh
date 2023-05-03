#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 20                                   # number of cores
#SBATCH --mem 150000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

#OrthoFinder v2.5.4
source package fc91613f-1095-4f67-b5aa-b86d702b36da

out_dir=../scratch/orthofinder

#Make output directory
mkdir -p ${out_dir}

orthofinder 	-f ../data/ncbi_data/proteins/gaeumannomyces/ \
		-t ${SLURM_CPUS_PER_TASK} \
		-o ${out_dir}/gaeumannomyces/

#Make directory for HOGs
mkdir -p ${out_dir}/gaeumannomyces/Results_Apr28/Single_Copy_HOGs

#Create FASTA sequence files for HOGs present in all taxa
singularity exec ~/programmes/orthofinder/230428_orthofinder.img create_files_for_hogs.sh \
	${out_dir}/gaeumannomyces/Results_Apr28/ ${out_dir}/gaeumannomyces/Results_Apr28/Single_Copy_HOGs/ N0

#Filter for single copy HOGs
for hog in $(ls ${out_dir}/gaeumannomyces/Results_Apr28/Single_Copy_HOGs/N0/HOG_Sequences/N0.HOG*)
do
        if [[ $(grep '>' $hog | sed 's/_EIv1.*//' | sed 's/\.faa_.*//' | uniq | wc -l) != 11 ]] || [[ $(grep '>' $hog | uniq | wc -l) != 11 ]]
        then
                rm $hog
        fi
done
