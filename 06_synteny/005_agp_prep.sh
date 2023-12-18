#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 1GB	                        # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

strain=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
dir=$(ls ../results/hifiasm_assemblies/*/v1.1/sequences/*${strain}* | sed -n 1p | sed 's|\(.*\)/.*|\1|')

awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${dir}/*${strain}_EI_v1.1_p_ctg.fa.fai > ${dir}/${strain}_EI_v1.1_p_ctg.tsv

awk '{print $1 "\t" $2 "\t" $3 "\t" 1 "\t" "W" "\t" $1 "\t" $2 "\t" $3 "\t" "+"}' ${dir}/${strain}_EI_v1.1_p_ctg.tsv > tmp${strain} && mv tmp${strain} ${dir}/${strain}_EI_v1.1_p_ctg.tsv

contig=($(awk -v var=${strain} '$1 == var' pseudochromosomes.tsv | awk '{print $2}'))
chr=($(awk -v var=${strain} '$1 == var' pseudochromosomes.tsv | awk '{print $3}' | sed 's/pseu_//'))
len=$(echo ${#contig[@]})

for (( i=0; i<$len; i++ ))
do

        sed -i "s/${contig[i]}/${chr[i]}/" ${dir}/${strain}_EI_v1.1_p_ctg.tsv

done

sed -i '/^chr/! s/[^ptg]*_ptg/ptg/' ${dir}/${strain}_EI_v1.1_p_ctg.tsv
sed -i 's/\t0\t/\t1\t/g' ${dir}/${strain}_EI_v1.1_p_ctg.tsv
