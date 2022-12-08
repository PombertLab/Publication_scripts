#!/usr/bin/bash
## On Spartacus
chown -R tomcat:jpombert /media/Data_1/apollo/
chown -R jpombert:jpombert /media/Data_1/apollo/data
mkdir /media/Data_1/apollo/data
mkdir /media/Data_1/apollo/data/E_hellem_50604

## Loading Fasta
/home/jpombert/Downloads/Apollo-2.4.1/web-app/jbrowse/bin/prepare-refseqs.pl \
--fasta /media/Data_4/jpombert/Microsporidia/E_hellem_50604/WebApollo/E_hellem_basecalled_v3_FINAL.fasta \
--out /media/Data_1/apollo/data/E_hellem_50604

# http://216.47.151.222:8085/apollo/annotator/index
# Organism -> Add new organism
# Name -> E_hellem_50604; Genus -> Encephalitozoon; Species -> hellem;
# Directory -> /media/Data_1/apollo/data/E_hellem_50604

/home/jpombert/Downloads/Apollo-2.4.1/docs/web_services/examples/groovy/add_organism.groovy \
-name E_hellem_50604 \
-url http://localhost:8085/apollo/ \
-directory /media/Data_1/apollo/data/E_hellem_50604 \
-username jpombert@iit.edu \
-password 'xxx'

## Prodigal
prodigal -c -f gff -i ../E_hellem_basecalled_v3_FINAL.fasta -o E_hellem_basecalled_v3_FINAL.gff
../../../scripts/splitGFF3.pl E_hellem_basecalled_v3_FINAL.gff

for k in {01..11}; do /home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl \
--gff /media/Data_4/jpombert/Microsporidia/E_hellem_50604/WebApollo/prodigal/chromosome_$k.gff3 \
--type CDS --subfeatureClasses '{"CDS": "orange-80pct"}' \
--trackLabel PRODIGAL \
--out /media/Data_1/apollo/data/E_hellem_50604; done

## TBLASTN
makeblastdb -in ../../E_hellem_basecalled_v3_FINAL.fasta -dbtype nucl -out Ehel_50604
tblastn -num_threads 10 -query Eint_50506_GCF_000146465.1_ASM14646v1_protein.faa -db DB/Ehel_50604 -outfmt 6 -out Eint_50506_GCF_000146465.1_ASM14646v1_protein.tblastn.6
tblastn -num_threads 10 -query Ehel_50504_GCF_000277815.2_ASM27781v3_protein.faa -db DB/Ehel_50604 -outfmt 6 -out Ehel_50504_GCF_000277815.2_ASM27781v3_protein.tblastn.6
tblastn -num_threads 10 -query Erom_SJ2008_GCF_000280035.1_ASM28003v2_protein.faa -db DB/Ehel_50604  -outfmt 6 -out Erom_SJ2008_GCF_000280035.1_ASM28003v2_protein.tblastn.6

## Loading TBLASTN
../../../scripts/getProducts.pl *.faa
../../../scripts/TBLASTN_to_GFF3.pl *.tblastn.6
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_hellem_50604/WebApollo/BLAST/Eint_50506_GCF_000146465.1_ASM14646v1_protein.gff \
--type match,match_part --subfeatureClasses '{"match_part": "orange-80pct"}' --trackLabel Eint_50506_tblastn --out /media/Data_1/apollo/data/E_hellem_50604
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_hellem_50604/WebApollo/BLAST/Ehel_50504_GCF_000277815.2_ASM27781v3_protein.gff \
--type match,match_part --subfeatureClasses '{"match_part": "green-80pct"}' --trackLabel Ehel_50504_tblastn --out /media/Data_1/apollo/data/E_hellem_50604
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_hellem_50604/WebApollo/BLAST/Erom_SJ2008_GCF_000280035.1_ASM28003v2_protein.gff \
--type match,match_part --subfeatureClasses '{"match_part": "blue-80pct"}' --trackLabel Erom_SJ2008_tblastn --out /media/Data_1/apollo/data/E_hellem_50604

## tRNAscan-2.0
../../../scripts/run_tRNAscan.pl ../E_hellem_basecalled_v3_FINAL.fasta
../../../scripts/tRNAscan_to_GFF3.pl *.tRNAs
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl \
--gff /media/Data_4/jpombert/Microsporidia/E_hellem_50604/WebApollo/tRNAscan/E_hellem_basecalled_v3_FINAL.fasta.tRNAs.gff \
--type tRNA \
--trackLabel tRNAscan-SE \
--out /media/Data_1/apollo/data/E_hellem_50604

## Loading BLASTN
blastn -num_threads 10 -query Ehel_50504_GCF_000277815.2_ASM27781v3_rna.fna \
-db DB/Ehel_50604 -outfmt 6 -out Ehel_50504_GCF_000277815.2_ASM27781v3_rna.blastn.6
blastn -num_threads 10 -query Eint_50506_GCF_000146465.1_ASM14646v1_rna.fna \
-db DB/Ehel_50604 -outfmt 6 -out Eint_50506_GCF_000146465.1_ASM14646v1_rna.blastn.6
../../../scripts/getProducts.pl *.fna
../../../scripts/TBLASTN_to_GFF3.pl *.blastn.6

/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl \
--gff /media/Data_4/jpombert/Microsporidia/E_hellem_50604/WebApollo/BLAST/Ehel_50504_GCF_000277815.2_ASM27781v3_rna.gff \
--type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' \
--trackLabel Ehel_50504_blastn --out /media/Data_1/apollo/data/E_hellem_50604

/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl \
--gff /media/Data_4/jpombert/Microsporidia/E_hellem_50604/WebApollo/BLAST/Eint_50506_GCF_000146465.1_ASM14646v1_rna.gff \
--type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' \
--trackLabel Eint_50506_blastn --out /media/Data_1/apollo/data/E_hellem_50604

## Creating BLAT Databases
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod a+x faToTwoBit
mkdir /media/Data_1/apollo/data/E_hellem_50604/twoBit
faToTwoBit /media/Data_4/jpombert/Microsporidia/E_hellem_50604/WebApollo/E_hellem_basecalled_v3_FINAL.fasta /media/Data_1/apollo/data/E_hellem_50604/twoBit/Ehel_50604.2bit
## In apollo
Blat database = /media/Data_1/apollo/data/E_hellem_50604/twoBit/Ehel_50604.2bit
