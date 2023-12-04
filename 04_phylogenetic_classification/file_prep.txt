#!/bin/sh

#Example one-liners to prepare fasta sequence headers in gene alignments to be identical prior to concatenation

#Removes accession and everything after key identifiers (e.g. name and voucher) based on certain words
sed 's/>[^ ]* />/' gene_alignment.fasta | sed 's/ glyce.*//' | sed 's/ GAPDH.*//' > tmp && mv tmp gene_alignment.fasta

#Replaces spaces, colons and brackets with underscores
sed 's/ /_/g' gene_alignment.fasta | sed 's/_$//' | sed 's/:/_/g' | sed 's/(//g' | sed 's/)//g' > tmp && mv tmp gene_alignment.fasta
