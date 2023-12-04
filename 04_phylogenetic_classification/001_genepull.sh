#!/bin/bash
#SBATCH -p ei-short                             # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 1                                    # number of cores
#SBATCH --mem 100                               # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=end,fail                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

source bedtools-2.30.0
source blast-2.7.1

out_dir=extract_genes

for strain in Gt-19d1 Gt-23d Gt-3aA1 Gt-4e Gt-8d Gt-CB1 Gt14LH10
do
		
	assembly=$(ls ../scratch/hifiasm-assemblies/${strain}/${strain}.asm.bp.p_ctg.fa)

	for marker in $(cat markers)
	do

		#If that marker is being used for the lineage...
		if [[ -f "${out_dir}/${marker}_example.fasta" ]]
	        then

			#Pull all hits for each marker from the assemblies
			echo "a" | ~/scripts/GenePull -g ${out_dir}/${marker}_example.fasta -a ${assembly} -o ${out_dir}/${strain}_${marker}

			#Fix headers
			sed -i "/^>/ s|>.*|>${strain}|" ${out_dir}/${strain}_${marker}.fa

			#If the marker is multi-copy...
			if [[ $(grep ">" ${out_dir}/${strain}_${marker}.fa | wc -l) > 1 ]]
			then
					
				#Filter out short hits
				length=$(grep ${marker} extract_genes/gene_lengths | awk '{print $2}')
				awk -v n=${length} '/^>/{ if(l>n) print b; b=$0;l=0;next } {l+=length;b=b ORS $0}END{if(l>n) print b }' ${out_dir}/${strain}_${marker}.fa > tmp && mv tmp ${out_dir}/${strain}_${marker}.fa

				#Make headers unique
                                awk 'BEGIN{RS=">";OFS="\n"}(NR>1){print ">"$1"_"(NR-1)"\n";$1="";print $0}' ${out_dir}/${strain}_${marker}.fa | awk '$0' > tmp && mv tmp ${out_dir}/${strain}_${marker}.fa
		
			fi

		fi

	done

done
