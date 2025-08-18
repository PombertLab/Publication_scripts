#!/usr/bin/env python3

name = 'round_robin_annotate.py'
version = '0.1.0'
updated = '2025-08-12'

from sys import argv
from subprocess import run
from shutil import which
from os.path import exists, isfile, basename
from os import makedirs
from argparse import ArgumentParser

############################################################
## Readme
############################################################

usage = f"""
NAME        {name}
VERSION     {version}
UPDATED     {updated}
SYNOPSIS    Assigns annotations to a set of query protein sequences based on homology identifed via a BLASTP search
            to a reference protein set.

USAGE       {name} \\
              -f 50507/50507.proteins \\
              -r REF/50506.faa

OPTIONS
-q (--query)        File containing protein sequences to be annotated
-f (--fastas)       File containing reference protein sequences (annotations must be in header, i.e. >E_int_50507 hexokinase)
-o (--outdir)       Output directory [Default:'RoundRobin']
"""

if len(argv) < 2:
    print(usage)
    exit()

############################################################
## Core functions
############################################################

def check_for_program(program:str=None) -> bool:

    """
    Checks for the presence of required program within the PATH; stops if missing.
    """

    if not which(program):

        raise Exception(f"{program} could not be found on PATH!")


def BLAST(query_file:str, sub_file:str, outdir:str) -> str:

    """
    Runs DIAMOND homology searches between the query and subject sequences 
    """

    query_basename = basename(query_file).split("/")[-1].split(".")[0]
    subject_basename = basename(sub_file).split("/")[-1].split(".")[0]
    outfile = f"{outdir}/{query_basename}_vs_{subject_basename}.diamond.blastp.6"

    print(f"BLASTing {query_basename} against {subject_basename}")

    if not isfile(query_file):
        raise Exception(f"Query file '{query_file}' does not exist.")
    if not isfile(sub_file):
        raise Exception(f"Subject file '{sub_file}' does not exist.")
    
    process = run(['diamond','blastp',
        '--query',f"{query_file}",
        '--db',f"{sub_file}",
        '--max-hsps','1',
        '--max-target-seqs','1',
        '--outfmt','6',
        '--out',f"{outfile}"
    ],capture_output=True)

    return outfile


def parse_protein_file(prot:str=None) -> tuple[set,dict]:

    """
    Parses FASTA protein files and grab loci and thier description from
    the FASTA headers [lines starting with >]. 
    """

    loci = set()
    annotations = {}

    with open(prot,'r') as IN:
        for line in IN:
            line = line.strip()
            if line[0] != '>':
                continue
            data = line[1:].split()
            loci.add(data[0])
            annotations[data[0]] = " ".join(data[1:]).split("[")[0].strip()

    return loci, annotations


def parse_BLAST_file(blast:str=None,invert_results:bool=False):

    """
    Parses the output of BLAST/DIAMOND outfmt 6 files; returns the corresponding matches
    """

    hits = {}

    with open(blast,'r') as IN:
        for line in IN:
            line = line.strip()
            query,subject = line.split("\t")[0:2]
            if invert_results:
                hits[subject] = query
            else:
                hits[query] = subject

    return hits


def check_for_bidirectionality(annotations:dict=None,sub_annotations:dict=None,label:str=None,loci:set=None,fh:dict=None,rh:dict=None) -> dict:

    """
    Checks if BLAST/DIAMOND homology hits are uni or bidirectional
    """

    annotations[label] = {}

    fhk = set(fh.keys())
    rhk = set(rh.keys())

    for locus in loci:

        direct = None
        annot = None
        match_loci = None

        if locus not in fhk and locus not in rhk:
            direct = "=/="
            annot = 'hypothetical protein'
            match_loci = 'N/A'
        elif locus in fhk and locus not in rhk:
            direct = "==>"
            annot = sub_annotations[fh[locus]]
            match_loci = fh[locus]
        elif locus not in fhk and locus in rhk:
            direct = "<=="
            annot = sub_annotations[rh[locus]]
            match_loci = rh[locus]
        else:
            direct = "<=>"
            annot = sub_annotations[fh[locus]]
            match_loci = fh[locus]

        annotations[label][locus] = {'direct':direct,'annot':annot,'locus':match_loci}

    return annotations


def annotate(query_file:str=None,loci:set=None,fastas:list=None,outdir:str=None):

    """
    Annotates bidirectional hits found by BLAST/DIAMOND homolgy searches.
    """

    annotations = {}

    for sub_file in fastas:

        subject_basename = basename(sub_file).split("/")[-1].split(".")[0]
        annotations[subject_basename] = {}

        _,sub_annotations = parse_protein_file(prot=sub_file)

        forward_hits = parse_BLAST_file(BLAST(query_file=query_file,sub_file=sub_file,outdir=outdir))
        
        reverse_hits = parse_BLAST_file(BLAST(sub_file=query_file,query_file=sub_file,outdir=outdir),invert_results=True)
        
        annotations = check_for_bidirectionality(
                        annotations=annotations,
                        sub_annotations=sub_annotations,
                        label=subject_basename,
                        loci=loci,
                        fh=forward_hits,
                        rh=reverse_hits
                    )

    return annotations


def write_out_annotations(loci:set=None,query_basename:str=None,annotations:dict=None,introns:dict=None,outdir:str=None):

    """
    Writes results to file ending with the .annotations extension.
    """

    MASTER_OUT = open(f"{outdir}/{query_basename}.annotations",'w')
    intron_keys = set([x for x in introns.keys()])

    try:
        for locus in sorted(loci):
            MASTER_OUT.write(f"{locus}")
            for label in annotations.keys():
                annotation = annotations[label][locus]['annot']
                directionality = annotations[label][locus]['direct']
                hit = annotations[label][locus]['locus']
                MASTER_OUT.write(f"\t{annotation}\t{directionality}\t{hit}")
                if hit in intron_keys:
                    MASTER_OUT.write(f"\t{introns[hit]}")
            MASTER_OUT.write("\n")
    finally:
        MASTER_OUT.close()

    return


def identify_introns(gff_file:str=None) -> dict:

    """
    Searches for the presence of introns in proteins by counting the number of exons
    in the corresponding CDS entries.
    """

    introns = {}

    with open(gff_file,'r') as IN:
        
        intron = None

        for line in IN:
            line = line.strip()
            
            if line == "":
                continue

            if line[0:11] == "/protein_id":
                introns[line.split()[-1].split("=")[-1].replace('"',"")] = intron
                continue

            if line[0:3] == "CDS":
                intron = line.count("..")
                intron = intron-1 if intron > 1 else None
                continue

    return introns


def main(args:ArgumentParser.parse_known_args=None) -> int:

    """
    Runs the script.
    """

    query_file=args.query
    fastas=args.fastas
    outdir=args.outdir

    makedirs(outdir,exist_ok=True)
    query_basename = basename(query_file).split(".")[0]

    loci,_ = parse_protein_file(prot=query_file)

    annotations = {}
    sub_annotations = {}

    annotations = annotate(query_file=query_file,loci=loci,fastas=fastas,outdir=outdir)

    introns = identify_introns(args.introns)

    write_out_annotations(loci=loci,query_basename=query_basename,annotations=annotations,introns=introns,outdir=outdir)

    return 0

############################################################
## Main
############################################################

if __name__ == '__main__':

    GetOptions = ArgumentParser()

    GetOptions.add_argument('-q','--query',required=True,type=str)
    GetOptions.add_argument('-f','--fastas',nargs='+',required=True,type=str)
    GetOptions.add_argument('-o','--outdir',default='RoundRobin',type=str)

    main(args=GetOptions.parse_known_args()[0])
