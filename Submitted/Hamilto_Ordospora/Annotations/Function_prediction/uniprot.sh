#!/usr/bin/bash

## Generating lists of all proteins queries; useful to identify proteins without hits in homology searches
./get_queries.pl *.fasta

## Running searches against Swiss-PROT
echo 'Querying BEOM2 against SwissProt...'
echo 'BEOM2 SwissProt searches started on:' >> swissprot.log; date >> swissprot.log
blastp -num_threads 10 -query BEOM2.proteins.fasta -db /media/Data_2/Uniprot/DB/SPROT -evalue 1e-10 -culling_limit 1 -outfmt 6 -out BEOM2.sprot.blastp.6
echo 'BEOM2 SwissProt searches completed on:' >> swissprot.log; date >> swissprot.log
./parse_UniProt_BLASTs.pl BEOM2.proteins.queries /media/Data_2/Uniprot/uniprot_sprot.products.hash BEOM2.sprot.blastp.6

echo 'Querying FIOER33 against SwissProt...'
echo 'FIOER33 SwissProt searches started on:' >> swissprot.log; date >> swissprot.log
blastp -num_threads 10 -query FIOER33.proteins.fasta -db /media/Data_2/Uniprot/DB/SPROT -evalue 1e-10 -culling_limit 1 -outfmt 6 -out FIOER33.sprot.blastp.6
echo 'FIOER33 SwissProt searches completed on:' >> swissprot.log; date >> swissprot.log
./parse_UniProt_BLASTs.pl FIOER33.proteins.queries /media/Data_2/Uniprot/uniprot_sprot.products.hash FIOER33.sprot.blastp.6

echo 'Querying ILG3 against SwissProt...'
echo 'ILG3 SwissProt searches started on:' >> swissprot.log; date >> swissprot.log
blastp -num_threads 10 -query ILG3.proteins.fasta -db /media/Data_2/Uniprot/DB/SPROT -evalue 1e-10 -culling_limit 1 -outfmt 6 -out ILG3.sprot.blastp.6
echo 'ILG3 SwissProt searches completed on:' >> swissprot.log; date >> swissprot.log
./parse_UniProt_BLASTs.pl ILG3.proteins.queries /media/Data_2/Uniprot/uniprot_sprot.products.hash ILG3.sprot.blastp.6

echo 'Querying ILBN2 against SwissProt...'
echo 'ILBN2 SwissProt searches started on:' >> swissprot.log; date >> swissprot.log
blastp -num_threads 10 -query ILBN2.proteins.fasta -db /media/Data_2/Uniprot/DB/SPROT -evalue 1e-10 -culling_limit 1 -outfmt 6 -out ILBN2.sprot.blastp.6
echo 'ILBN2 SwissProt searches completed on:' >> swissprot.log; date >> swissprot.log
./parse_UniProt_BLASTs.pl ILBN2.proteins.queries /media/Data_2/Uniprot/uniprot_sprot.products.hash ILBN2.sprot.blastp.6

## Running searches against TREMBL
echo 'Querying BEOM2 against TREMBL...'
echo 'BEOM2 TREMBL searches started on:' >> trembl.log; date >> trembl.log
blastp -num_threads 10 -query BEOM2.proteins.fasta -db /media/Data_2/Uniprot/DB/TREMBL -evalue 1e-10 -culling_limit 1 -outfmt 6 -out BEOM2.trembl.blastp.6
echo 'BEOM2 TREMBL searches completed on:' >> trembl.log; date >> trembl.log
./parse_UniProt_BLASTs.pl BEOM2.proteins.queries /media/Data_2/Uniprot/uniprot_trembl.products.hash BEOM2.trembl.blastp.6

echo 'Querying FIOER33 against TREMBL...'
echo 'FIOER33 SwissProt searches started on:' >> trembl.log; date >> trembl.log
blastp -num_threads 10 -query FIOER33.proteins.fasta -db /media/Data_2/Uniprot/DB/TREMBL -evalue 1e-10 -culling_limit 1 -outfmt 6 -out FIOER33.trembl.blastp.6
echo 'FIOER33 TREMBL searches completed on:' >> trembl.log; date >> trembl.log
./parse_UniProt_BLASTs.pl FIOER33.proteins.queries /media/Data_2/Uniprot/uniprot_trembl.products.hash FIOER33.trembl.blastp.6

echo 'Querying ILG3 against TREMBL...'
echo 'ILG3 TREMBL searches started on:' >> trembl.log; date >> trembl.log
blastp -num_threads 10 -query ILG3.proteins.fasta -db /media/Data_2/Uniprot/DB/TREMBL -evalue 1e-10 -culling_limit 1 -outfmt 6 -out ILG3.trembl.blastp.6
echo 'ILG3 TREMBL searches completed on:' >> trembl.log; date >> trembl.log
./parse_UniProt_BLASTs.pl ILG3.proteins.queries /media/Data_2/Uniprot/uniprot_trembl.products.hash ILG3.trembl.blastp.6

echo 'Querying ILBN2 against TREMBL...'
echo 'ILBN2 TREMBL searches started on:' >> trembl.log; date >> trembl.log
blastp -num_threads 10 -query ILBN2.proteins.fasta -db /media/Data_2/Uniprot/DB/TREMBL -evalue 1e-10 -culling_limit 1 -outfmt 6 -out ILBN2.trembl.blastp.6
echo 'ILBN2 TREMBL searches completed on:' >> trembl.log; date >> trembl.log
./parse_UniProt_BLASTs.pl ILBN2.proteins.queries /media/Data_2/Uniprot/uniprot_trembl.products.hash ILBN2.trembl.blastp.6

