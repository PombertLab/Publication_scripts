#!/usr/bin/bash

./get_GBM1_NCBI_annotations \
-g GCF_000091225.1_ASM9122v1_genomic.gff \
-o GBM1_NCBI_annotations.tsv

./get_GBM1_MSDB_annotations \
-g MicrosporidiaDB-56_EcuniculiGBM1.gff \
-o GBM1_MSDB_annotations.tsv

$TDFI/Prediction/AlphaFold2/parse_af_results.pl \
-a 3DFI/Folding/ALPHAFOLD_3D_Parsed/ \
-o AlphaFold_pLDDTs.tsv

mkdir VOROCNN
mkdir VOROCNN/ALPHAFOLD/
cd VOROCNN/ALPHAFOLD

vorocnn -i 3DFI/Folding/ALPHAFOLD_3D_Parsed -v $VORONATA

cd ../
mkdir VOROCNN/RAPTORX
cd VOROCNN/RAPTORX

vorocnn -i 3DFI/Folding/RAPTORX/PDB -v $VORONATA

cd ../

../vorocnn_average.pl \
-i RAPTORX/vorocnn_output/*.scores ALPHAFOLD/vorocnn_output/*.scores \
-o vorocnn_averages.tsv

$TDFI/Homology_search/run_MICAN_on_GESAMT_results.pl \
-t 3DFI \
-r $TDFI_DB/RCSB_PDB $TDFI_DB/RCSB_PDB_obsolete \
-o MICAN

./create_master_table.pl \
-n GBM1_NCBI_annotations.tsv \
-d GBM1_MSDB_annotations.tsv \
-p AlphaFold_pLDDTs.tsv \
-g 3DFI/Homology/GESAMT/All_GESAMT_matches_per_protein.tsv \
-m MICAN/MICAN_raw.tsv \
-v VOROCNN/vorocnn_averages.tsv > MASTER_TABLE.tsv

./get_score_distribution.py plddt MASTER_TABLE.tsv

./get_score_distribution.py q MASTER_TABLE.tsv

./get_score_distribution.py tm MASTER_TABLE.tsv