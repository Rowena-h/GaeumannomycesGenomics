#!/bin/bash

cd ../scratch/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Analysis

for f in {Gt-19d1,Gt-4e,Gh-2C17,Gh-1B17,Gt-3aA1,Gt-8d,Gt-CB1,Gt14LH10,Gt-23d}
do 
	mkdir -p ${f}/Analysis/eggnog-mapper-2.1.9_CBG
	cd  ${f}/Analysis/eggnog-mapper-2.1.9_CBG
	ln -s /ei/.project-scratch/d/d2c0bfb1-c37a-4211-9493-86b15d4e773e/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/${f}/${f}_EIv1.release.gff3.pep.fasta ./
	sbatch -p ei-long -J ${f}-eggnog --mem=40G -c 16 -o eggnog_%j_%N.out -e eggnog_%j_%N.err  --mail-type=end,fail --mail-user=olivera@nbia.ac.uk --wrap "source eggnog-mapper-2.1.9_CBG && /usr/bin/time -v emapper.py --cpu 8 --itype proteins -m diamond --cpu 8  --pfam_realign none --go_evidence non-electronic --report_orthologs -d /ei/cb/common/Databases/eggnog/eggnog_proteins.dmnd --report_no_hits  --override -i ${f}_EIv1.release.gff3.pep.fasta -o ${f}_eggnog"
mark 
done
