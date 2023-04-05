#!/bin/bash
#SBATCH -p ei-medium                      	# queue
#SBATCH -N 1                            	# number of nodes
#SBATCH -c 10                            	# number of cores
#SBATCH --mem 400000                     	# memory pool for all cores
#SBATCH --output ../scratch/logs/%x.o%j         # STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=rowena.hill@earlham.ac.uk	# send-to address

source minimap2-2.21

strain_file=$(awk '{print $2}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)

#Convert unitig GFA file to fasta
awk '/^S/{print ">"$2;print $3}' ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.bp.p_utg.gfa > ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.bp.p_utg.fa

read_file=$(awk '{print $3}' ../strains | sed -n ${SLURM_ARRAY_TASK_ID}p)
readsmap_dir=../scratch/007_readsmap/${strain_file}

#Map reads to unitigs
minimap2 \
        -ax map-hifi -t ${SLURM_CPUS_PER_TASK} \
        ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_utg.fa \
        ../scratch/hifi-reads/${strain_file}.fastq.gz > ${readsmap_dir}/${strain_file}_utg_aln.sam

source samtools-1.13

#Convert to indexed, sorted BAM file
samtools view -S -b ${readsmap_dir}/${STRAIN_FILE}_utg_aln.sam > ${readsmap_dir}/${STRAIN_FILE}_utg_aln.bam
samtools sort ${readsmap_dir}/${STRAIN_FILE}_utg_aln.bam -o ${readsmap_dir}/${STRAIN_FILE}_utg_aln_sorted.bam
samtools index ${readsmap_dir}/${STRAIN_FILE}_utg_aln_sorted.bam

source blast-2.10

blastn_dir=../scratch/006_blastn/${strain_file}

#BLAST unitigs against nt database
blastn \
    -task megablast \
    -db /ei/public/databases/blast/ncbi/nt_20210521/nt \
    -query ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.bp.p_utg.fa \
    -out ${blastn_dir}/${strain_file}_utg_hits \
    -outfmt "6 qseqid staxids bitscore std" \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -evalue 1e-25

source blobtools-1.0.1

blobtools_dir=../scratch/008_blobtools/${strain_file}

for rank in species family order
do

        blobtools create -i ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.bp.p_utg.fa \
			 -b ${readsmap_dir}/${strain_file}_utg_aln_sorted.bam \
			 -t ${blastn_dir}/${strain_file}_utg_hits \
			 -o ${blobtools_dir}/${strain_file}_utg

        blobtools view -i ${blobtools_dir}/${strain_file}_utg.blobDB.json -o ${blobtools_dir}/${rank} -r ${rank}

        blobtools plot -i ${blobtools_dir}/${strain_file}_utg.blobDB.json -o ${blobtools_dir}/${rank} -r ${rank}

done

source seqtk-1.0

out_dir=../scratch/009_remove_contam/${strain_file}

#Make output directory
mkdir -p ${out_dir}

if [[ ${strain_file} = "PG3bAQG1a" ]]
then

        #Split mixed sample
        awk 'NR==FNR {A[$1]++;next} $6 in A {print $1}' PG3bAQG1a_paraphaeosphaeria_families ${blobtools_dir}/family.${strain_file}_utg.blobDB.table.txt | \
        seqtk subseq ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_utg.fa - > ${out_dir}_para/${strain_file}_para.asm.bp.p_utg.fa

        awk 'NR==FNR {A[$1]++;next} $6 in A {print $1}' PG3bAQG1a_pyrenophora_families ${blobtools_dir}/family.${strain_file}_utg.blobDB.table.txt | \
        seqtk subseq ../scratch/hifiasm-assemblies/${strain_file}/${strain_file}.asm.bp.p_utg.fa - > ${out_dir}_pyreno/${strain_file}_pyreno.asm.bp.p_utg.fa

fi
