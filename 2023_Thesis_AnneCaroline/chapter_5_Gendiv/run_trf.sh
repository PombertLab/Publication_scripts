#!/usr/bin/env bash

OUTDIR="./Repeats/TRF"

## outdir
mkdir -p $OUTDIR

## trf
trf \
  Hellem_T2T_ATCC_50451.fna \
  2 \
  5 \
  7 \
  80 \
  10 \
  50 \
  2000 \
  1

# ./trf409.linux64 File Match Mismatch Delta PM PI Minscore MaxPeriod [options]
# File = Hellem_T2T_ATCC_50451.fna
# match weight 2;
# mismatch penalty 5;
# indel penalty 7,
# match probability 80;
# indel probability 10;
# min score 50;
# max period 2000;
# max tr length 1

## cleanup
mv *.html $OUTDIR/