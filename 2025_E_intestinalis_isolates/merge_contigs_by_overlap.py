#!/usr/bin/env python3

name = 'merge_contigs_by_overlap.py'
version = '0.2.0'
updated = '2025-08-12'

from sys import argv
from argparse import ArgumentParser
from os import makedirs
from os.path import basename,isdir
from textwrap import wrap

############################################################
## Readme
############################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Attempts to merge contigs by searching for overlaps.

USAGE       {name} \\
			-f ORIENTED_ASSEMBLY/50507.oriented.fasta \\
			-m ORIENTED_ASSEMBLY/all.map \\
			-o MERGED_ASSEMBLY

OPTIONS
-f (--fasta)    FASTA file of contigs
-m (--map)      Chromosome map file
-p (--prefix)   Contig prefix to assign [Default: 'contig_']
-o (--outdir)   Output directory [Default: MERGED_CONTIGS]
"""

if len(argv) < 2:
    print(usage)
    exit()

############################################################
## Command lines arguments
############################################################

GetOptions = ArgumentParser()

GetOptions.add_argument("-f","--fasta",required=True)
GetOptions.add_argument("-m","--map",required=True)
GetOptions.add_argument("-p","--prefix",default='contig_')
GetOptions.add_argument("-o","--outdir",default="MERGED_CONTIGS")

args = GetOptions.parse_args()

fasta = args.fasta
mapp = args.map
prefix = args.prefix
outdir = args.outdir

############################################################
## Temp dir
############################################################

filename = basename(fasta).split(".")[0]
temp_dir = f"{outdir}/{filename}"

if not isdir(temp_dir):
    makedirs(temp_dir,mode=0o755)

############################################################
## Loading FASTA file
############################################################

sequences = {}
locus = False

FASTA = open(fasta,'r')
for line in FASTA:
    
    line = line.strip()

    if line[0] == '>':
        
        locus = line[1:]
        sequences[locus] = ""
    
    elif locus:

        sequences[locus] += line

FASTA.close()

############################################################
## Read contig chromosome mappings
############################################################

IN = open(mapp,'r')

locus = ""
hits = {}

for line in IN:

    line = line.strip()
    
    ## Skip empty and header lines
    if line == "" or line[0:2] == "##":
        continue
    
    ## Get chromosome accessions
    if line[0:2] == ">>":
        locus = line[2:].split("\t")[0]
        hits[locus] = {}
        continue

    contig,match_type,_,qstart,qend,_,_,sstart,send = line[1:].split("\t")[:9]

    ## Skip secondary mappings
    if match_type != 'Primary':
        continue

    hits[locus][contig] = {'qstart':int(qstart),'qend':int(qend),'sstart':int(sstart),'send':int(send)}

IN.close()

############################################################
## Read contig chromosome mappings
############################################################

OVERLAP = open(f"{temp_dir}/overlap.alignment",'w')
LOG = open(f"{temp_dir}/merge.log",'w')
LOG.write(f"## > CONTIG\tOVERLAPS\tMATCHED\tOVERLAP_BPS\tPIDENT\tACTION\n")

merged_sequences = {}
merged_to = {}

## Iterate over
for ref in hits.keys():

    LOG.write(f"Reference= {ref}\n\n")

    ## References with 1 mapping contig will have no overlaps
    if len(hits[ref].keys()) < 2:
        key = "".join([x for x in hits[ref].keys()])
        merged_sequences[key] = sequences[key]
        LOG.write(f"> {key}\n\n")
        continue

    previous_start = False
    previous_qstart = False
    previous_end = False
    previous_contig = False

    for sub in sorted(hits[ref].keys(),key = lambda x: hits[ref][x]['sstart']):

        start = hits[ref][sub]['sstart']
        end = hits[ref][sub]['send']

        ## First mapping contig will not have an overlap on 5' end
        if not previous_start:

            merged_sequences[sub] = sequences[sub]
            LOG.write(f"> {sub}\t-\t-\t-\t-\t-\n")

            previous_qstart,previous_start,previous_end,previous_contig = [qstart,start,end,sub]

            continue

        ## Current contig starts mapping after previous contig ends
        if start > previous_end:

            if previous_contig not in merged_sequences.keys() and previous_contig not in merged_to.keys():
                merged_sequences[previous_contig] = sequences[previous_contig]
                LOG.write(f"> {previous_contig}\t-\t-\t-\t-\t-\n")

            if sub not in merged_sequences.keys():
                merged_sequences[sub] = sequences[sub]
                LOG.write(f"> {sub}\t-\t-\t-\t-\t-\n")

            previous_qstart,previous_start,previous_end,previous_contig = [qstart,start,end,sub]

            continue

        ## Current contig ends before previous mapping contig
        if end < previous_end:

            LOG.write(f"> {sub}\t completely encompassed by {previous_contig}n")
            continue

        overlap = previous_end - start + 1
        
        seq1 = sequences[previous_contig][len(sequences[previous_contig])-overlap:]
        seq2 = sequences[sub][:overlap]
        match = []

        for index,x in enumerate(seq1):
            if x == seq2[index]:
                match.append("|")
            else:
                match.append(".")

        match = "".join(match)

        matched = sum([1 for x in match if x == "|"]) 
        pident = matched/overlap*100

        merged_contig = previous_contig

        ## Gets the first contig where mapping overlap occurs
        while merged_contig in merged_to.keys():
            
            merged_contig = merged_to[merged_contig]

        
        if pident == 100:

            merged_sequences[merged_contig] += sequences[sub][overlap:]
            merged_to[sub] = previous_contig

            LOG.write(f"> {sub}\t{merged_contig}\t{overlap}\t{matched}\t{pident:.2f}%\tMERGED\n")

        else:
            print("Overlap found, but not high enough bases matched to merge!")
            merged_sequences[previous_contig] = sequences[previous_contig]
            LOG.write(f"> {sub}\t{merged_contig}\t{overlap}\t{matched}\t{pident:.2f}%\tREJECTED_MERGE\n")
        

        seq1 = wrap(seq1,60)
        seq2 = wrap(seq2,60)
        match = wrap(match,60)

        buffer = len(str(previous_end))

        OVERLAP.write(f"Query= {previous_contig}\n")
        OVERLAP.write(f"> {sub}\n\n")
        for x in range(len(seq1)):
            OVERLAP.write(f"Query  {previous_end-overlap+(60*x)+1:<{buffer}d}   {seq1[x]} {previous_end-overlap+(60*x)+len(seq1[x]):<{buffer}d}\n")
            temp_match = match[x].replace("."," ")
            OVERLAP.write(f"       {' '*buffer}   {temp_match}\n")
            OVERLAP.write(f"Sbjct  {(60*x)+1:<{buffer}d}   {seq2[x]} {(60*x)+len(seq1[x]):<{buffer}d}\n\n")
        
        OVERLAP.write(f"\n")

        previous_qstart,previous_start,previous_end,previous_contig = [qstart,start,end,sub]
    
    LOG.write(f"\n")

OVERLAP.close()

############################################################
## Write merged contigs
############################################################

OUT = open(f"{temp_dir}/{filename}.merged.fasta",'w')

buffer = len(str(len(merged_sequences.keys())))
contig_count = 1

for contig in sorted(merged_sequences.keys(),key = lambda x: len(merged_sequences[x]),reverse=True):

    OUT.write(f">{contig}\n")

    seq = "\n".join(wrap(merged_sequences[contig],60))

    OUT.write(f"{seq}\n")

    contig_count+= 1

OUT.close()