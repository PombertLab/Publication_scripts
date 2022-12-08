#!/usr/bin/bash

ODIR=READ_DISTRIBUTIONS

COLOR=darkorange
XLIM=120000
TICKS=5
HEIGHT=4

# cuniculi
DATA=E_cuniculi_50602.nanopore
~/GitHub/Misc/read_len_plot.py \
  -f $DATA.fastq.gz \
  -d $ODIR \
  -o $DATA.svg $DATA.png $DATA.pdf \
  -c $COLOR \
  -h $HEIGHT \
  -x $XLIM \
  -t $TICKS

# hellem
DATA=E_hellem_50604.nanopore
~/GitHub/Misc/read_len_plot.py \
  -f $DATA.fastq.gz \
  -d $ODIR \
  -o $DATA.svg $DATA.png $DATA.pdf \
  -c $COLOR \
  -h $HEIGHT \
  -x $XLIM \
  -t $TICKS


# intestinalis
DATA=E_intestinalis_50506.nanopore
~/GitHub/Misc/read_len_plot.py \
  -f $DATA.fastq.gz \
  -d $ODIR \
  -o $DATA.svg $DATA.png $DATA.pdf \
  -c $COLOR \
  -h $HEIGHT \
  -x $XLIM \
  -t $TICKS


DATA=E_intestinalis_50506.pacbio
~/GitHub/Misc/read_len_plot.py \
  -f $DATA.fastq.gz \
  -d $ODIR \
  -o $DATA.svg $DATA.png $DATA.pdf \
  -c $COLOR \
  -h $HEIGHT \
  -x $XLIM \
  -t $TICKS
