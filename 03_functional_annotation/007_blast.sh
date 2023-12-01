#!/bin/bash
#SBATCH -p ei-medium				# queue
#SBATCH -N 1					# number of nodes
#SBATCH -c 1					# number of cores
#SBATCH --mem 1G				# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

strain=$(awk '{print $1}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
proteins=../data/ncbi_data/proteins/gaeumannomyces/${strain}_EIv1.release.gff3.pep.repisoform.fasta
out_dir=../scratch/functional_annotation

#BLAST CSEPs against PHI-base

mkdir -p ${out_dir}/004_CSEP_prediction

source seqtk-1.0

seqtk subseq ${proteins} ${strain}_CSEPs.txt > ${out_dir}/004_CSEP_prediction/${strain}_CSEPs.faa

source blast-2.10

blastp \
	-query ${out_dir}/004_CSEP_prediction/${strain}_CSEPs.faa \
	-db phi-base_current.fas \
	-outfmt "6 qseqid sseqid evalue bitscore pident length" \
	-evalue 1e-25 \
	-out ${out_dir}/004_CSEP_prediction/${strain}_phibase_blast.tsv \
	-num_threads ${SLURM_CPUS_PER_TASK}

#Filter for top bitscore result per gene
sort -r -n -k4 < ${out_dir}/004_CSEP_prediction/${strain}_phibase_blast.tsv | awk '!x[$1]++' ${out_dir}/004_CSEP_prediction/${strain}_phibase_blast.tsv | sort -k 1b,1  > ${out_dir}/004_CSEP_prediction/${strain}_phibase_blast_tophits.tsv
#Add data to CSEPs file
join -a1 -a2 -t $'\t' -o 1.1 2.2 2.3 2.4 2.5 -1 1 -2 1 ${strain}_CSEPs.txt ${out_dir}/004_CSEP_prediction/${strain}_phibase_blast_tophits.tsv > ${strain}_CSEPs.tsv
awk -F "\t" '{gsub(/\#/,"\t",$2);print $0}' ${strain}_CSEPs.tsv | sed 's/ /\t/g' > ${strain}tmp && mv ${strain}tmp ${strain}_CSEPs.tsv

#BLAST for avenacinase gene and MAT loci

mkdir -p ${out_dir}/005_avenacinase

makeblastdb \
	-in ${proteins} \
	-dbtype prot

blastx \
	-query avenacinase_reference.fa \
	-db ${proteins}	\
	-outfmt "6 qseqid sseqid evalue bitscore pident length" \
	-evalue 1e-25 \
	-out ${out_dir}/005_avenacinase/${strain}_avenacinase.tsv \
	-num_threads ${SLURM_CPUS_PER_TASK}

mkdir -p ${out_dir}/005_MAT

blastp \
	-query MAT_reference.fa \
	-db ${proteins} \
	-outfmt "6 qseqid sseqid evalue bitscore pident sstart send length" \
	-evalue 1e-25 \
	-out ${out_dir}/005_MAT/${strain}_MAT.tsv \
	-num_threads ${SLURM_CPUS_PER_TASK}
