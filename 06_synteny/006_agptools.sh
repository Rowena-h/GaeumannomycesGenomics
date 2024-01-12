#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 1GB	                        # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

strain=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
organism=$(awk '{print $3}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
new_strain=$(awk '{print $4}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
dir=$(ls ../results/hifiasm_assemblies/*/v1.1/sequences/*${strain}* | sed -n 1p | sed 's|\(.*\)/.*|\1|')
out_dir=$(../results/ncbi_submission/${strain})

#Make output directory
mkdir -p ${out_dir}


## ASSEMBLY ##

#Merge into chromosomes
singularity exec ~/programmes/agptools/agptools.img agptools assemble \
	${dir}/*${strain}_EI_v1.1_p_ctg.fa \
	${out_dir}/${strain}_EI_v1.1_p_ctg_curated.agp > ${out_dir}/${strain}_merged.fa

#Rename headers
for chr in $(grep ">chr" ${out_dir}/${strain}_merged.fa | sed 's/>//')
do
	num=$(echo $chr | sed 's/chr//')	
	sed -i "s/${chr}/${chr} [organism=${organism}] [strain=${new_strain}] [location=chromosome] [chromosome=${num}]/" ${out_dir}/${strain}_merged.fa
done	

sed -i "/>ptg/ s/$/ [organism=${organism}] [strain=${new_strain}]/" ${out_dir}/${strain}_merged.fa


## ANNOTATION ##

awk '/^[^#]/ {print $0}' ../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/*/${strain}_EIv1.release.gff3 > tmp${strain}.gff

#Flip feature orientations for inverted contigs
contigs=$(awk '{print $1}' tmp${strain}.gff | sort | uniq)
flip=$(awk '$9 == "-" {print $6}' ${out_dir}/${strain}_EI_v1.1_p_ctg_curated.agp)

rm tmp${strain}inverted.gff

for contig in $contigs
do
	if [[ $flip == *"$contig"* ]]
	then
		awk -v var="$contig" 'BEGIN{ FS = OFS = "\t" } $1 == var{sub("-", "plus", $7); sub("+", "minus", $7); print $0}' tmp${strain}.gff >> tmp${strain}inverted.gff
	else
		awk -v var="$contig" 'BEGIN{ FS = OFS = "\t" } $1 == var{print $0}' tmp${strain}.gff >> tmp${strain}inverted.gff
	fi
done

awk 'BEGIN{ FS = OFS = "\t" } {sub("plus", "+", $7); sub("minus", "-", $7); print $0}' tmp${strain}inverted.gff \
	> tmp${strain}inverted2.gff && mv tmp${strain}inverted2.gff tmp${strain}inverted.gff

awk '/^[^#]/ {print $1 "\t" $4 "\t" $5}' tmp${strain}inverted.gff > tmp${strain}inverted.bed

#Transform feature coordinates to match merged assembly
~/programmes/agptools/agptools.img agptools transform \
	tmp${strain}inverted.bed \
	${out_dir}/${strain}_EI_v1.1_p_ctg_curated.agp \
	-o tmp${strain}invertedmerged.bed

#Add new coordinates to GFF
awk 'BEGIN{ FS = OFS = "\t" } FNR==NR{a[NR]=$1;b[NR]=$2;c[NR]=$3;next}{$1=a[FNR];$4=b[FNR];$5=c[FNR]}1' \
	tmp${strain}invertedmerged.bed tmp${strain}inverted.gff > ${out_dir}/${strain}_merged.gff3

rm tmp${strain}invertedmerged.bed tmp${strain}inverted.gff tmp${strain}inverted.bed tmp${strain}.gff
