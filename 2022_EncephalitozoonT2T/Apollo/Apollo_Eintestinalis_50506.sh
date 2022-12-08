#!/usr/bin/bash
## On MiniMe; /media/Data_1/jpombert/Microsporidia/E_intestinalis_50506/WebApollo

## On Spartacus
chown -R tomcat:jpombert /media/Data_1/apollo/
chown -R jpombert:jpombert /media/Data_1/apollo/data
mkdir /media/Data_1/apollo/data
mkdir /media/Data_1apollo/data/E_intestinalis_50506

## Loading Fasta
/home/jpombert/Downloads/Apollo-2.4.1/web-app/jbrowse/bin/prepare-refseqs.pl \
--fasta /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/Eintestinalis_basecalled_v2.fasta \
--out /media/Data_1/apollo/data/E_intestinalis_50506

# http://216.47.151.222:8085/apollo/annotator/index
# Organism -> Add new organism
# Name -> E_intestinalis_50506; Genus -> Encephalitozoon; Species -> intestinalis; Directory -> /media/Data_1/apollo/data/E_intestinalis_50506

## Prodigal
prodigal -c -f gff -i ../Eintestinalis_basecalled_v2.fasta -o Eintestinalis_basecalled_v2.gff
../scripts/splitGFF3.pl Eintestinalis_basecalled_v2.gff

for k in {01..11}; do /home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl \
--gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/prodigal/chromosome_$k.gff3 \
--type CDS --subfeatureClasses '{"CDS": "orange-80pct"}' \
--trackLabel PRODIGAL \
--out /media/Data_1/apollo/data/E_intestinalis_50506; done

## BAM file illumina_PE Minimap2; BAM URL is bam folder inside /media/Data_1/apollo/data/E_intestinalis_50506, doesn't work with full path
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/add-bam-track.pl \
--bam_url bam/Eintestinalis_basecalled_v2.fasta.minimap2.bam \ 
--label illumina_PE \
--key "illumina_PE" \
-i /media/Data_1/apollo/data/E_intestinalis_50506/trackList.json

## TBLASTN
makeblastdb -in Eintestinalis_basecalled_v2.fasta -dbtype nucl -out Eint_50506
tblastn -num_threads 10 -query Eint_50506_GCF_000146465.1_ASM14646v1_protein.faa -db DB/Eint_50506 -outfmt 6 -out Eint_50506_GCF_000146465.1_ASM14646v1_protein.tblastn.6
tblastn -num_threads 10 -query Ehel_50504_GCF_000277815.2_ASM27781v3_protein.faa -db DB/Eint_50506 -outfmt 6 -out Ehel_50504_GCF_000277815.2_ASM27781v3_protein.tblastn.6
tblastn -num_threads 10 -query Erom_SJ2008_GCF_000280035.1_ASM28003v2_protein.faa -db DB/Eint_50506  -outfmt 6 -out Erom_SJ2008_GCF_000280035.1_ASM28003v2_protein.tblastn.6

## Loading TBLASTN
./getProducts.pl *.faa
./TBLASTN_to_GFF3.pl *.tblastn.6
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/BLAST/Eint_50506_GCF_000146465.1_ASM14646v1_protein.gff --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' --trackLabel Eint_50506_tblastn --out /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/BLAST/Ehel_50504_GCF_000277815.2_ASM27781v3_protein.gff --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' --trackLabel Ehel_50504_tblastn --out /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/BLAST/Erom_SJ2008_GCF_000280035.1_ASM28003v2_protein.gff --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' --trackLabel Erom_SJ2008_tblastn --out /media/Data_1/apollo/data/E_intestinalis_50506

## tRNAscan-2.0
./scripts/run_tRNAscan.pl Eintestinalis_basecalled_v2.fasta
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/tRNAscan/Eintestinalis_basecalled_v2.fasta.tRNAs.gff --type tRNA --trackLabel tRNAscan-SE --out /media/Data_1/apollo/data/E_intestinalis_50506

## RNAmmer doesn't seem to work on microsporidia!

## Loading BLASTN
blastn -num_threads 10 -query Eint_50506_GCF_000146465.1_ASM14646v1_rna.fna -db DB/Eint_50506 -outfmt 6 -out Eint_50506_GCF_000146465.1_ASM14646v1_rna.blastn.6
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/BLAST/Eint_50506_GCF_000146465.1_ASM14646v1_rna.gff --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' --trackLabel Eint_50506_blastn --out /media/Data_1/apollo/data/E_intestinalis_50506

## Creating BLAT Databases
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod a+x faToTwoBit
mkdir /media/Data_1/apollo/data/E_intestinalis_50506/twoBit
faToTwoBit /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/Eintestinalis_basecalled_v2.fasta /media/Data_1/apollo/data/E_intestinalis_50506/twoBit/Eint_50506.2bit
## In apollo
Blat database = /media/Data_1/apollo/data/E_intestinalis_50506/twoBit/Eint_50506.2bit

## Reverse-complement
## Don't forget to reverse complement chromosomes to get the same order as other Encephalitozoon spp.!
## Removing tracks, then  redo with chromosomes in the right order for easier comparisons!
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/remove-track.pl --delete --trackLabel Eint_50506_tblastn --dir /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/remove-track.pl --delete --trackLabel Ehel_50504_tblastn --dir /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/remove-track.pl --delete --trackLabel Erom_SJ2008_tblastn --dir /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/remove-track.pl --delete --trackLabel tRNAscan-SE --dir /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/remove-track.pl --delete --trackLabel PRODIGAL --dir /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/remove-track.pl --delete --trackLabel illumina_PE --dir /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/remove-track.pl --delete --trackLabel Eint_50506_blastn --dir /media/Data_1/apollo/data/E_intestinalis_50506

## Reverse-complementing with reseq from EMBOSS; then redo with the right sequences
splitMultifasta.pl ../Eintestinalis_basecalled_v2.fasta
for k in {01,03,04,05,06,08,09,11}; do revseq chromosome_$k.fsa chromosome_${k}_rev.fsa; done
mkdir WRONG_ORIENTATION
for k in {01,03,04,05,06,08,09,11}; do mv chromosome_$k.fsa WRONG_ORIENTATION/; done
cat *.fsa > E_intestinalis_50506.fasta

## Loading Fasta
/home/jpombert/Downloads/Apollo-2.4.1/web-app/jbrowse/bin/prepare-refseqs.pl \
--fasta /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/E_intestinalis_50506.fasta \
--out /media/Data_1/apollo/data/E_intestinalis_50506

faToTwoBit /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/E_intestinalis_50506.fasta /media/Data_1/apollo/data/E_intestinalis_50506/twoBit/Eint_50506.2bit

## Prodigal
prodigal -c -f gff -i ../E_intestinalis_50506.fasta -o E_intestinalis_50506.gff
../scripts/splitGFF3.pl E_intestinalis_50506.gff

for k in {01..11}; do /home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl \
--gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/prodigal/chromosome_$k.gff3 \
--type CDS --subfeatureClasses '{"CDS": "orange-80pct"}' \
--trackLabel PRODIGAL \
--out /media/Data_1/apollo/data/E_intestinalis_50506; done

## Loading BLASTN
blastn -num_threads 10 -query Eint_50506_GCF_000146465.1_ASM14646v1_rna.fna -db DB/Eint_50506 -outfmt 6 -out Eint_50506_GCF_000146465.1_ASM14646v1_rna.blastn.6
../scripts/getProducts.pl *.fna
../scripts/TBLASTN_to_GFF3.pl *.blastn.6
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/BLAST/Eint_50506_GCF_000146465.1_ASM14646v1_rna.gff --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' --trackLabel Eint_50506_blastn --out /media/Data_1/apollo/data/E_intestinalis_50506

## TBLASTN
makeblastdb -in E_intestinalis_50506.fasta -dbtype nucl -out Eint_50506
tblastn -num_threads 10 -query Eint_50506_GCF_000146465.1_ASM14646v1_protein.faa -db DB/Eint_50506 -outfmt 6 -out Eint_50506_GCF_000146465.1_ASM14646v1_protein.tblastn.6
tblastn -num_threads 10 -query Ehel_50504_GCF_000277815.2_ASM27781v3_protein.faa -db DB/Eint_50506 -outfmt 6 -out Ehel_50504_GCF_000277815.2_ASM27781v3_protein.tblastn.6
tblastn -num_threads 10 -query Erom_SJ2008_GCF_000280035.1_ASM28003v2_protein.faa -db DB/Eint_50506  -outfmt 6 -out Erom_SJ2008_GCF_000280035.1_ASM28003v2_protein.tblastn.6

## Loading TBLASTN
../scripts/getProducts.pl *.faa
../scripts/TBLASTN_to_GFF3.pl *.tblastn.6
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/BLAST/Eint_50506_GCF_000146465.1_ASM14646v1_protein.gff --type match,match_part --subfeatureClasses '{"match_part": "orange-80pct"}' --trackLabel Eint_50506_tblastn --out /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/BLAST/Ehel_50504_GCF_000277815.2_ASM27781v3_protein.gff --type match,match_part --subfeatureClasses '{"match_part": "green-80pct"}' --trackLabel Ehel_50504_tblastn --out /media/Data_1/apollo/data/E_intestinalis_50506
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/BLAST/Erom_SJ2008_GCF_000280035.1_ASM28003v2_protein.gff --type match,match_part --subfeatureClasses '{"match_part": "blue-80pct"}' --trackLabel Erom_SJ2008_tblastn --out /media/Data_1/apollo/data/E_intestinalis_50506

## tRNAscan-2.0
../scripts/run_tRNAscan.pl E_intestinalis_50506.fasta
../scripts/scripts/tRNAscan_to_GFF3.pl *.tRNAs
/home/jpombert/Downloads/Apollo-2.4.1/jbrowse/bin/flatfile-to-json.pl --gff /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/tRNAscan/E_intestinalis_50506.fasta.tRNAs.gff --type tRNA --trackLabel tRNAscan-SE --out /media/Data_1/apollo/data/E_intestinalis_50506

## Creating BLAT Databases
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod a+x faToTwoBit
mkdir /media/Data_1/apollo/data/E_intestinalis_50506/twoBit
faToTwoBit /media/Data_4/jpombert/Microsporidia/E_intestinalis_50506/WebApollo/E_intestinalis_50506.fasta /media/Data_1/apollo/data/E_intestinalis_50506/twoBit/Eint_50506.2bit
## In apollo
Blat database = /media/Data_1/apollo/data/E_intestinalis_50506/twoBit/Eint_50506.2bit
