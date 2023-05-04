#!/usr/bin/bash

## Performed on nanopore data basecalled with Guppy 6.0+

./MinIONQC.R \
  -i sequencing_summary_E_cuniculi_50602.txt \
  -o E_cuniculi_50602 \
  -f 'png'

./MinIONQC.R \
  -i sequencing_summary_E_hellem_50604.txt \
  -o E_hellem_50604 \
  -f 'png'


./MinIONQC.R \
  -i sequencing_summary_E_intestinalis_50506.txt \
  -o E_intestinalis_50506 \
  -f 'png'

