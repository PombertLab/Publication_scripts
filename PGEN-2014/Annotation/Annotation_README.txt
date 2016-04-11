Running maker
********************************************************************
1.	Sort contigs so that the minimum size is 500 bp:
	sort_contigs.pl -m 500 fasta_input_file fasta_output_file
	sort_contigs.pl -z -b input_500bp+.fasta input_500bp+_sorted.fasta	## sort the contigs from longuest to smallest

2.	Rename the contigs to that they are named by +1 increments
	number_contigs.pl input_500bp+.fasta

3.	Run Maker on the sorted file

********************************************************************

PFAMs searches
********************************************************************
1.	catenate all proteins into a multifasta file
2.	perl headers_new.pl file.fasta	## Rename with proper BioProject ID
3.	do a batch PFAM search (http://pfam.sanger.ac.uk/search#tabview=tab1)	## may need to split in more than one query
4.	download hitdata.txt file
5.	perl hitdata_evalue_filter.pl hitdata.txt ## filter according to desired evalue
6.	May need to remove whitespaces with whitespaces_to_tab.pl
7.	Generate a list containing locus_tag, PFAM, evalue, then add proper NCBI description tag (must be manually curated)
********************************************************************

Annotating
********************************************************************
1.	verify that all contigs are in one .gff file (some are in two files gff.ann and gff.seq)
2.	GFF3_to_FASTA.pl *.gff
3.	fasta_to_string.pl *.fsa
4.	extractGFF_features_v3.0.pl *.gff
5.	perl ../Scripts/verify_completedness_v2.pl *.feat
6.	perl ../Scripts/features_to_TBL_v8.0.pl *.feat
7.	perl ../Scripts/partial_mRNAs_tags_to_tbl.pl *.tbl
8.	perl ../Scripts/partial_CDS_tags_to_tbl_v3.pl *.ends
9.	perl ../Scripts/Add_annotated_products.pl *.tbl3	## Based on PFAM searches (evalue threshold used = 1e-30)
10.	rm *.tbl	## remove tmp tbl files
11.	rm *.tbl2	## remove tmp tbl files
12.	rm *.tbl3	## remove tmp tbl files
13.	perl ../Scripts/correctLocus.pl *.tbl4
14.	tbl2asn -t template.sbt -p TBL2ASN/ -g
********************************************************************