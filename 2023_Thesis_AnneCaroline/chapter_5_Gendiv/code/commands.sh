#!/usr/bin/env bash

## 
./coding_percentage.pl -f Hellem_T2T_ATCC_50451.gbff -o Circos -c CDS

##
./nucleotide_biases.pl -f Hellem_T2T_ATCC_50451.fna -o Circos -c Biases

## Kmers with circos output option
./k_counter.py -i Hellem_T2T_ATCC_50451.fna -o kmer_redundancy -t 12 -k 7 -s 1000 -w 500 

## BLAST
FASTA=Hellem_T2T_ATCC_50451.fna
SHORT=50451
blastn -query $FASTA -subject $FASTA -outfmt 6 -evalue 1e-10 -out $SHORT.blastn.6
./b2links.pl -b $SHORT.blastn.6  -c

circos -conf circos.conf