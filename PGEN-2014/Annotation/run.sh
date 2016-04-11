#!/usr/bin/bash
cd GFFs/
perl ../Scripts/verify_completedness_v2.pl *.feat
perl ../Scripts/features_to_TBL_v8.0.pl *.feat
perl ../Scripts/partial_mRNAs_tags_to_tbl.pl *.tbl
perl ../Scripts/partial_CDS_tags_to_tbl_v3.pl *.ends
perl ../Scripts/Add_annotated_products.pl *.tbl3	## Based on PFAM searches (evalue threshold used = 1e-30)
rm *.tbl	## remove tmp tbl files
rm *.tbl2	## remove tmp tbl files
rm *.tbl3	## remove tmp tbl files
perl ../Scripts/correctLocus.pl *.tbl4
rm *.tbl4
mv *.tbl ../FILES/
cd ../
tbl2asn -p FILES/ -t template.sbt -w assembly.cmt -M n -Z discrep