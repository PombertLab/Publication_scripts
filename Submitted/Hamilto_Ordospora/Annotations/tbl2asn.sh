#!/usr/bin/bash

echo "Converting EMBL to TBL ..."
cd GFF3/
/media/Data_1/jpombert/WebApollo_Karen/scripts/EMBLtoTBL.pl *.embl
cd ../
echo "Running TBL2ASN..."
tbl2asn -t template.sbt -w genome.asm -p GFF3/ -g -M n -Z discrep -V b -H 12/31/2018
echo "Concatenating validations .val files into val.txt..."
cat GFF3/*.val > val.txt

## To generate a single SQN:
# tbl2asn -t template.sbt -w genome.asm -p GFF3/ -g -M n -Z discrep -H 12/31/2018 -o combined_contigs_BEOM2.sqn


