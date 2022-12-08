#!/usr/bin/bash

# busco --list-datasets


# eukaryota_odb10
# - fungi_odb10
#          - ascomycota_odb10
#              - dothideomycetes_odb10
#                  - capnodiales_odb10
#                  - pleosporales_odb10
#              - eurotiomycetes_odb10
#                  - chaetothyriales_odb10
#                  - eurotiales_odb10
#                  - onygenales_odb10
#              - leotiomycetes_odb10
#                  - helotiales_odb10
#              - saccharomycetes_odb10
#              - sordariomycetes_odb10
#                  - glomerellales_odb10
#                  - hypocreales_odb10
#          - basidiomycota_odb10
#              - agaricomycetes_odb10
#                  - agaricales_odb10
#                  - boletales_odb10
#                  - polyporales_odb10
#              - tremellomycetes_odb10
#          - microsporidia_odb10
#          - mucoromycota_odb10
#              - mucorales_odb10

## busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]

# E_cuniculi_50602
ENC=E_cuniculi_50602
busco \
  -i $ENC.prot \
  -l fungi_odb10 \
  -o ${ENC}_vs_Fungi \
  -m proteins

busco \
  -i $ENC.prot \
  -l microsporidia_odb10 \
  -o ${ENC}_vs_Microsporidia \
  -m proteins

busco \
  -i $ENC.prot \
  -l eukaryota_odb10 \
  -o ${ENC}_vs_Eukaryota \
  -m proteins

# E_hellem_50604
ENC=E_hellem_50604
busco \
  -i $ENC.prot \
  -l fungi_odb10 \
  -o ${ENC}_vs_Fungi \
  -m proteins

busco \
  -i $ENC.prot \
  -l microsporidia_odb10 \
  -o ${ENC}_vs_Microsporidia \
  -m proteins

busco \
  -i $ENC.prot \
  -l eukaryota_odb10 \
  -o ${ENC}_vs_Eukaryota \
  -m proteins

# E_intestinalis_50506
ENC=E_intestinalis_50506
busco \
  -i $ENC.prot \
  -l fungi_odb10 \
  -o ${ENC}_vs_Fungi \
  -m proteins

busco \
  -i $ENC.prot \
  -l microsporidia_odb10 \
  -o ${ENC}_vs_Microsporidia \
  -m proteins

busco \
  -i $ENC.prot \
  -l eukaryota_odb10 \
  -o ${ENC}_vs_Eukaryota \
  -m proteins
