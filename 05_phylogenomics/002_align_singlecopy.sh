#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 1000                              # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

orthofinder_dir=$(ls -td ../scratch/orthofinder/gaeumannomyces/* | head -1)
#num_orthos=$(ls ${orthofinder_dir}/Single_Copy_HOGs/N0/HOG_Sequences/*.fa | wc -l)

out_dir=alignments/orthogroups

#Make output directory
mkdir -p ${out_dir}

#Total array required = 7066
#values=({1..5000})
#values=({5001..7066})

#Submit array job to align each single copy HOGs
sbatch -p ei-medium -N 1 -c 1 --mem 150000 -a 0-2065 --output ../scratch/logs/gaeumannomyces_align.o%j --mail-type=FAIL --mail-user=rowena.hill@earlham.ac.uk --wrap="source mafft-7.271; values=({5001..7066}); array_value=\${values[\$SLURM_ARRAY_TASK_ID]}; ortho=\$(ls \"${orthofinder_dir}\"/Single_Copy_HOGs/N0/HOG_Sequences/*.fa | sed 's#^.*/##' | sed 's#\.fa##' | sed -n \${array_value}p); mafft \"${orthofinder_dir}\"/Single_Copy_HOGs/N0/HOG_Sequences/\${ortho}.fa > alignments/orthogroups/\${ortho}_aln.fa; sed -i '/>/ s/\(.*\)_.*$/\1/' alignments/orthogroups/\${ortho}_aln.fa; singularity exec ~/programmes/trimAl/trimAl.img trimal -in alignments/orthogroups/\${ortho}_aln.fa -fasta -gappyout > alignments/orthogroups/\${ortho}_aln_trimmed.fa"
