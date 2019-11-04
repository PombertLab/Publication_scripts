### Converting GFF3 to EMBL format ###
## Splitting the WebApollo GFF3 file per contig
gunzip BEOM2.gff3.gz
splitGFF3.pl -g BEOM2.gff3

## Converting WebApollo GFF3 annotations to EBML format
WebApolloGFF3toEMBL.pl -p CWI36 -g contig*.gff3

## Manual checking with Artemis to ensure that features were transferred
art contig-0013.embl

## Exporting protein from EMBL files; useful if proteins are added/removed after WebApollo curation, otherwise redundant with WebApolloGFF3toEMBL.pl
EMBLtoPROT *.embl

## Replacing locus_tag prefixes if temporary ones were used while waiting for the NCBI ones
replace_prefixes.pl NEWPREFIX *.embl



### Function prediction ### 
## Predicting functions with InterProScan 5
bash interproscan.sh

## Downloading the SwissProt/UniProt databases
get_UniProt.pl -s -t -n 20 -l download.log

## Creating tab-delimited lists of sequence in the SwissProt/UniProt databases
hash_uniprot.pl uniprot_sprot.fasta uniprot_trembl.fasta

## Running BLAST searches against SwissProt/UniProt
bash uniprot.sh

## Generating a list of all proteins queried, inncluding those with no hits against SwissProt/UniProt
get_queries.pl *.fasta

## Parsing the result of InterProScan 5 and SwissProt/UniProt searches
parse_annotators.pl -q BEOM2.proteins.queries -sl sprot.list -sb BEOM2.sprot.blastp.6 sprot.list -tl trembl.list -tb BEOM2.trembl.blastp.6 -ip BEOM2.interpro.tsv

## Curating the annotations
curate_annotations.pl -r -i BEOM2.annotations




### TBL2ASN ###
## Adding taxonomic info to FASTA (.fsa) files; Useful for TBL2ASN
add_taxinfo.pl -sp 'Hamiltosporidium magnivora' -is BE-OM-2 -fa *.fsa

## Converting EMBL to TBL
EMBLtoTBL.pl *.embl

## Generate a template.sbt file per genome
https://submit.ncbi.nlm.nih.gov/genbank/template/submission/

## Create a structure commments file (genome.asm) per genome; e.g.
StructuredCommentPrefix	##Genome-Assembly-Data-START##
Assembly Method	Ray v. AUGUST-2017
Assembly Name	Version_1
Long Assembly Name	BEOM2 v1
Genome Coverage	25.87x
Sequencing Technology	Illumina
StructuredCommentSuffix ##Genome-Assembly-Data-END##

## Running TBL2ASN
bash tbl2asn.sh

## Fixing misc errors Artemis
check_errors_1.sh
check_errors_2.sh
check_errors_3.sh

## Running TBL2ASN again until all errors are fixed...
bash tbl2asn.sh
