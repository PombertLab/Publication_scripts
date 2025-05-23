## WebApollo example (BEOM2)
## Creating data directories
mkdir /media/Data_2/apollo/data/BEOM2_NEW
chown -R tomcat:jpombert /media/Data_2/apollo/
chown -R jpombert:jpombert /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/prepare-refseqs.pl --fasta /media/Data_1/jpombert/WebApollo_Karen/BEOM2/BEOM2.contigs.fasta --out /media/Data_2/apollo/data/BEOM2_NEW

## Creatig organism in WebApollo
http://localhost:8085/apollo/annotator/index
Organism -> Add new organism
Name -> BEOM2; Genus -> Hamiltosporidium; Species -> magnivora; Directory -> /media/Data_2/apollo/data/BEOM2_NEW

## Loading MAKER annotation tracks
for k in {0001..3550}; do /home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/MAKER_GFFs/contig-$k.augustus.gff --type match,match_part --subfeatureClasses '{"match_part": "orange-80pct"}' --trackLabel MAKER_AUG_ROZELLA --out /media/Data_2/apollo/data/BEOM2_NEW; done
for k in {0001..3550}; do /home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/MAKER_GFFs/contig-$k.genemark.gff --type match,match_part --subfeatureClasses '{"match_part": "orange-80pct"}' --trackLabel MAKER_GENEMARK --out /media/Data_2/apollo/data/BEOM2_NEW; done
for k in {0001..3550}; do /home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/MAKER_GFFs/contig-$k.repeats.gff --type match,match_part --subfeatureClasses '{"match_part": "magenta-80pct"}' --trackLabel REPEATS --out /media/Data_2/apollo/data/BEOM2_NEW; done

## Loading GeneMark-ES predictions as user annotations; there is a bug in some cases where a few annotations are not loaded, use Artemis to fix that?
## We are missing 8 annotations; 5031 instead of 5039, fixing it manually is possible
for k in {0001..3550}; do /home/jpombert/Downloads/Apollo-2.0.6/tools/data/add_transcripts_from_gff3_to_annotations.pl -U http://localhost:8085/apollo -u jpombert@iit.edu -p 'password' -i /media/Data_1/jpombert/WebApollo_Karen/BEOM2/MAKER_GFFs/contig-$k.genemark.gff -t match -e match_part -o "BEOM2_NEW" 1>>output.log 2>> error.log; done

## Loading PRODIGAL > 200 aa
for k in {0001..3550}; do /home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/PRODIGAL/SPLIT_200/contig-$k.gff --type CDS --subfeatureClasses '{"CDS": "orange-80pct"}' --trackLabel PRODIGAL_200 --out /media/Data_2/apollo/data/BEOM2_NEW; done

## Loading TBLASTN searches in WebApollo
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/Rozella.proteins.gff --trackLabel Rozella --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' --out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-34_EcuniculiGBM1_AnnotatedProteins.gff --trackLabel EcuniculiGBM1 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' --out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_AalgeraePRA109_AnnotatedProteins.gff --trackLabel AalgeraePRA109 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_AalgeraePRA339_AnnotatedProteins.gff --trackLabel AalgeraePRA339 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_EaedisUSNM41457_AnnotatedProteins.gff --trackLabel EaedisUSNM41457 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_MdaphniaeUGP3_AnnotatedProteins.gff --trackLabel MdaphniaeUGP3 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_NausubeliERTm2_AnnotatedProteins.gff --trackLabel NausubeliERTm2 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_NausubeliERTm6_AnnotatedProteins.gff --trackLabel NausubeliERTm6 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_NbombycisCQ1_AnnotatedProteins.gff --trackLabel NbombycisCQ1 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_NceranaeBRL01_AnnotatedProteins.gff --trackLabel NceranaeBRL01 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_NparisiiERTm1_AnnotatedProteins.gff --trackLabel NparisiiERTm1 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_NparisiiERTm3_AnnotatedProteins.gff --trackLabel NparisiiERTm3 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_OcolligataOC4_AnnotatedProteins.gff --trackLabel OcolligataOC4 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_PneurophiliaMK1_AnnotatedProteins.gff --trackLabel PneurophiliaMK1 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_Slophii42_110_AnnotatedProteins.gff --trackLabel Slophii42 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_ThominisUnknown_AnnotatedProteins.gff --trackLabel ThominisUnknown --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_VcorneaeATCC50505_AnnotatedProteins.gff --trackLabel VcorneaeATCC50505 --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW
/home/jpombert/Downloads/Apollo-2.0.6/bin/flatfile-to-json.pl --gff /media/Data_1/jpombert/WebApollo_Karen/BEOM2/TBLASTN/MicrosporidiaDB-33_Vculicisfloridensis_AnnotatedProteins.gff --trackLabel Vculicisfloridensis --type match,match_part --subfeatureClasses '{"match_part": "brightgreen-80pct"}' -out /media/Data_2/apollo/data/BEOM2_NEW

## Perform manual curation of annotations
## Export curated annotations
Select 'Ref Sequence' tab
Export -> GFF3; select ALL; select GFF3 with FASTA; click Export
mv Annotations.gff3.gz BEOM2.gff3.gz
