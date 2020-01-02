#!/usr/bin/bash

## Running InterProScan 5 searches
echo 'BEOM2 InterProScan started on:' >> interproscan.log; date >> interproscan.log
/media/Data_2/opt/interproscan-5.26-65.0/interproscan.sh -cpu 10 -i BEOM2.proteins.fasta -iprlookup -goterms -pa -b BEOM2.interpro
echo 'BEOM2 InterProScan completed on:' >> interproscan.log; date >> interproscan.log

echo 'FIOER33 InterProScan started on:' >> interproscan.log; date >> interproscan.log
/media/Data_2/opt/interproscan-5.26-65.0/interproscan.sh -cpu 10 -i FIOER33.proteins.fasta -iprlookup -goterms -pa -b FIOER33.interpro
echo 'FIOER33 InterProScan completed on:' >> interproscan.log; date >> interproscan.log

echo 'ILG3 InterProScan started on:' >> interproscan.log; date >> interproscan.log
/media/Data_2/opt/interproscan-5.26-65.0/interproscan.sh -cpu 10 -i ILG3.proteins.fasta -iprlookup -goterms -pa -b ILG3.interpro
echo 'ILG3 InterProScan completed on:' >> interproscan.log; date >> interproscan.log

echo 'ILBN2 InterProScan started on:' >> interproscan.log; date >> interproscan.log
/media/Data_2/opt/interproscan-5.26-65.0/interproscan.sh -cpu 10 -i ILBN2.proteins.fasta -iprlookup -goterms -pa -b ILBN2.interpro
echo 'ILBN2 InterProScan completed on:' >> interproscan.log; date >> interproscan.log

