#!/usr/bin/env python3

name = 'process_fasta_sequences.py'
version = '0.2.2'
updated = '2025-08-12'

from sys import argv
from argparse import ArgumentParser
from textwrap import wrap
from os.path import isdir, basename
from os import makedirs

############################################################
## Readme
############################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Provides meaningful names to fasta sequences, orders them by size, and
            can remove sequences of undesirable length.

USAGE       {name} \\
              -f 50507.fasta \\
              -n 1000

OPTIONS
-f (--fasta)    Input fasta file(s)
-o (--outdir)   Output directory [Default: processed_fastas]
-p (--prefix)   Prefix to give contigs [Default: contig_]
-n (--min)      Minimum contig length [Default: 500]
-x (--max)      Maximum contig length [Default: None]
"""

if len(argv) < 2:
    print(usage)
    exit()

############################################################
## Command lines arguments
############################################################

GetOptions = ArgumentParser()

GetOptions.add_argument("-f","--fasta",nargs='+',required=True)
GetOptions.add_argument("-o","--outdir",default='processed_sequences')
GetOptions.add_argument("-p","--prefix",default='contig_')
GetOptions.add_argument("-n","--min",type=int,default=500)
GetOptions.add_argument("-x","--max",type=int,default=False)

args = GetOptions.parse_args()

fastas = args.fasta
outdir = args.outdir
prefix = args.prefix
lmin = args.min
lmax = args.max

############################################################
## Outdir
############################################################

if not isdir(outdir):
    makedirs(outdir,mode=0o755)

############################################################
## Iterating through FASTA files
############################################################

for file in fastas:

    print(f"\tProcessing {file}")

    sequences = {}
    locus = False
    filename = basename(file).split(".")[0]

    temp_dir = f"{outdir}/{filename}"

    if not isdir(temp_dir):
        makedirs(temp_dir,mode=0o755)

    ############################################################
    ## Read fasta file
    ############################################################

    IN = open(file,'r')

    for line in IN:
        
        line = line.strip()

        if line == "":
            continue

        if line[0] == '>':
            
            locus = line[1:]
            sequences[locus] = ""
        
        elif locus:

            sequences[locus] += line

    IN.close()

    ############################################################
    ## Size select sequences
    ############################################################

    keep = []

    for locus in sequences.keys():

        ## Discard sequences that are too short
        if len(sequences[locus]) <= lmin:
            continue

        ## Discard sequences that are too long (why, I don't know, but the customer is always right)
        if lmax and len(sequences[locus]) >= lmax:
            continue

        keep.append(locus)

    buffer = len(str(len(keep)))

    OUT = open(f"{temp_dir}/{filename}.processed.fasta",'w')
    LOG = open(f"{temp_dir}/contig_name_links.tsv",'w')
    LOG.write("## NEW_NAME\tOLD_NAME\n")
    seq_count = 1

    ############################################################
    ## Write out sequences
    ############################################################
    
    for key in sorted(keep,key=lambda x: len(sequences[x]),reverse=True):

        OUT.write(f">{prefix}{seq_count:0{buffer}d}\n")
        LOG.write(f"{prefix}{seq_count:0{buffer}d}\t{key}\n")
        seq = "\n".join(wrap(sequences[key],60))
        OUT.write(f"{seq}\n")
        seq_count += 1
    
    OUT.close()
    LOG.close()
