#!/bin/bash
#SBATCH -p ei-medium                            # queue
#SBATCH -N 1                                    # number of nodes
#SBATCH -c 2                                    # number of cores
#SBATCH --mem 150000                            # memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL                    # notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk   # send-to address

out_dir=starfish

#Make output directory
mkdir -p ${out_dir}


## Format files for starfish ##

mkdir -p ${out_dir}/inputs

#cp ../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/genome/*_EI_v1.1_p_ctg.fa inputs
#cp ../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/*/*_EIv1.release.gff3 inputs

#for file in $(ls ${out_dir}/inputs/*.fa)
#do
#	sed -i 's/_EI_v1\.1//g' $file
#done

#for file in $(ls ${out_dir}/inputs/*.gff3)
#do
#	sed -i 's/_EI_v1\.1//g' $file
#	sed -i '/^#/d' $file
#	awk '{if($3 == "gene") print $0}' $file > tmp && mv tmp $file
#done

realpath ${out_dir}/inputs/*.fa | perl -pe 's/^(.+?([^\/]+?).fa)$/\2\t\1/' > ${out_dir}/ome2assembly.txt
realpath ${out_dir}/inputs/*.gff3 | perl -pe 's/^(.+?([^\/]+?).gff3)$/\2\t\1/' > ${out_dir}/ome2gff.txt

cat $(awk '{print $2}' ${out_dir}/ome2gff.txt) > ${out_dir}/gaeumannomyces.gff3

mkdir -p ${out_dir}/blastdb

cut -f2 ${out_dir}/ome2assembly.txt | xargs cat > ${out_dir}/blastdb/gaeumannomyces.assemblies.fna

singularity exec ~/programmes/starfish/starfish_231101.img makeblastdb \
	-in ${out_dir}/blastdb/gaeumannomyces.assemblies.fna \
	-out ${out_dir}/blastdb/gaeumannomyces.assemblies \
	-parse_seqids \
	-dbtype nucl
	
cut -f1,5 $(awk '{print $2}' ../strains | \
	sed 's|^|\.\./scratch/functional_annotation/003_eggnog-mapper/|' | \
	sed 's|$|/*|') | \
	grep -v '#' | \
	perl -pe 's/(^.+?)\t.+,([^,]+)$/\1\t\2/' | \
	perl -pe 's/@/\t/' > ${out_dir}/gaeumannomyces.gene2og.txt

~/programmes/starfish/starfish/scripts/geneOG2mclFormat.pl \
	-i ${out_dir}/gaeumannomyces.gene2og.txt \
	-o ${out_dir}


## Run starfish ##

mkdir -p ${out_dir}/geneFinder

singularity exec ~/programmes/starfish/starfish_231101.img starfish annotate \
	-x gaeumannomyces_tyr \
	-a ${out_dir}/ome2assembly.txt \
	-g ${out_dir}/ome2gff.txt \
	-p ~/programmes/starfish/starfish/database/YRsuperfams.p1-512.hmm \
	-P ~/programmes/starfish/starfish/database/YRsuperfamRefs.faa \
	-i tyr \
	-o ${out_dir}/geneFinder/ \
	-T ${SLURM_CPUS_PER_TASK}

singularity exec ~/programmes/starfish/starfish_231101.img starfish consolidate \
	-o ${out_dir} \
	-g ${out_dir}/gaeumannomyces.gff3 \
	-G ${out_dir}/geneFinder/gaeumannomyces_tyr.filt_intersect.gff

realpath ${out_dir}/gaeumannomyces_tyr.filt_intersect.consolidated.gff | \
	perl -pe 's/^/gaeumannomyces\t/' > ${out_dir}/ome2consolidatedGFF.txt

singularity exec ~/programmes/starfish/starfish_231101.img starfish sketch \
	-m 10000 \
	-q ${out_dir}/geneFinder/gaeumannomyces_tyr.filt_intersect.ids \
	-g ${out_dir}/ome2consolidatedGFF.txt \
	-i s \
	-x gaeumannomyces \
	-o ${out_dir}/geneFinder/

grep -P '\ttyr\t' ${out_dir}/geneFinder/gaeumannomyces.bed > ${out_dir}/geneFinder/gaeumannomyces.tyr.bed 

mkdir -p ${out_dir}/elementFinder

singularity exec ~/programmes/starfish/starfish_231101.img starfish insert \
	-a ${out_dir}/ome2assembly.txt \
	-d ${out_dir}/blastdb/gaeumannomyces.assemblies \
	-b ${out_dir}/geneFinder/gaeumannomyces.tyr.bed \
	-i tyr \
	-x gaeumannomyces \
	-o ${out_dir}/elementFinder/ \
	-T ${SLURM_CPUS_PER_TASK}

singularity exec ~/programmes/starfish/starfish_231101.img starfish flank \
	-a ${out_dir}/ome2assembly.txt \
	-b ${out_dir}/elementFinder/gaeumannomyces.insert.bed \
	-x gaeumannomyces \
	-o ${out_dir}/elementFinder/

singularity exec ~/programmes/starfish/starfish_231101.img starfish summarize \
	-a ${out_dir}/ome2assembly.txt \
	-b ${out_dir}/elementFinder/gaeumannomyces.flank.bed \
	-x gaeumannomyces \
	-o ${out_dir}/elementFinder/ \
	-S ${out_dir}/elementFinder/gaeumannomyces.insert.stats \
	-f ${out_dir}/elementFinder/gaeumannomyces.flank.singleDR.stats \
	-g ${out_dir}/ome2consolidatedGFF.txt \
	-A ${out_dir}/gaeumannomyces.gene2og.txt \
	-t ${out_dir}/geneFinder/gaeumannomyces_tyr.filt_intersect.ids 

mkdir -p ${out_dir}/pairViz

singularity exec ~/programmes/starfish/starfish_231101.img starfish pair-viz \
	-m all \
	-t empty \
	-A nucmer \
	-a ${out_dir}/ome2assembly.txt \
	-b ${out_dir}/elementFinder/gaeumannomyces.elements.bed \
	-f ${out_dir}/elementFinder/gaeumannomyces.flank.singleDR.stats \
	-S ${out_dir}/elementFinder/gaeumannomyces.elements.named.stats \
	-o ${out_dir}/pairViz/ \
	-T ${SLURM_CPUS_PER_TASK}

mkdir -p ${out_dir}/regionFinder

singularity exec ~/programmes/starfish/starfish_231101.img mmseqs easy-cluster \
	${out_dir}/geneFinder/gaeumannomyces_tyr.filt_intersect.fas \
	${out_dir}/regionFinder/gaeumannomyces_tyr \
	${out_dir}/regionFinder/ \
	--min-seq-id 0.5 \
	-c 0.25 \
	--alignment-mode 3 \
	--cov-mode 0 \
	--cluster-reassign \
	--threads ${SLURM_CPUS_PER_TASK} \

~/programmes/starfish/starfish/scripts/mmseqs2mclFormat.pl \
	-i ${out_dir}/regionFinder/gaeumannomyces_tyr_cluster.tsv \
	-g fam \
	-o ${out_dir}/regionFinder/

singularity exec ~/programmes/starfish/starfish_231101.img starfish sim \
	-m element \
	-t nucl \
	-b ${out_dir}/elementFinder/gaeumannomyces.elements.bed \
	-x gaeumannomyces \
	-o ${out_dir}/regionFinder/ \
	-a ${out_dir}/ome2assembly.txt

singularity exec ~/programmes/starfish/starfish_231101.img starfish group \
	-m mcl \
	-s ${out_dir}/regionFinder/gaeumannomyces.element.nucl.sim \
	-i hap \
	-o ${out_dir}/regionFinder/ \
	-t 0.05

grep -P '\tcap\t' ${out_dir}/elementFinder/gaeumannomyces.elements.bed | cut -f4,7 > ${out_dir}/regionFinder/gaeumannomyces.cap2ship.txt

~/programmes/starfish/starfish/scripts/searchReplace.pl \
	-i ${out_dir}/regionFinder/gaeumannomyces_tyr_cluster.mcl \
	-r ${out_dir}/regionFinder/gaeumannomyces.cap2ship.txt > ${out_dir}/regionFinder/gaeumannomyces.element_cluster.mcl

~/programmes/starfish/starfish/scripts/mergeGroupfiles.pl \
	-t ${out_dir}/regionFinder/gaeumannomyces.element_cluster.mcl \
	-q ${out_dir}/regionFinder/gaeumannomyces.element.nucl.I1.5.mcl > ${out_dir}/regionFinder/gaeumannomyces.element.fam-hap.mcl

grep -f <(comm -23 <(cut -f1 ${out_dir}/geneFinder/gaeumannomyces_tyr.filt_intersect.ids | \
	sort) <(grep -P '\tcap\t|\ttyr\t' ${out_dir}/elementFinder/gaeumannomyces.elements.bed | \
	cut -f4| sort)) ${out_dir}/geneFinder/gaeumannomyces.tyr.bed > ${out_dir}/regionFinder/unaffiliated_tyrs.bed

~/programmes/starfish/starfish/scripts/filterOG.pl \
	-O ${out_dir}/gaeumannomyces.gene2og.mcl \
	-a 1 \
	-c 5 \
	-o ${out_dir}

singularity exec ~/programmes/starfish/starfish_231101.img starfish dereplicate \
	-e ${out_dir}/regionFinder/gaeumannomyces.element.fam-hap.mcl \
	-t ${out_dir}/regionFinder/unaffiliated_tyrs.bed \
	-F ${out_dir}/elementFinder/gaeumannomyces.elements.feat \
	-S ${out_dir}/elementFinder/gaeumannomyces.elements.named.stats \
	-O ${out_dir}/gaeumannomyces.gene2og.a1.c5.txt \
	-g ${out_dir}/ome2gff.txt \
	-x gaeumannomyces \
	-o ${out_dir}/regionFinder/ \
	--flanking 3