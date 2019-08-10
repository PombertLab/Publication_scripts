#!/usr/bin/perl

## Pombert Lab, IIT 2014
## Extracts open reading frames (ORFs) from fasta seqquences with Getorf from EMBOSS
## Requires the EMBOSS package

use strict;
use warnings;

## ORF type:
##	0	Translation of regions between STOP codons 
##	1	Translation of regions between START and STOP codons
##	2	Nucleic sequences between STOP codons
##	3	Nucleic sequences between START and STOP codons
##	4	Nucleotides flanking START codons
##	5	Nucleotides flanking initial STOP codons
##	6	Nucleotides flanking ending STOP codons

my $type = 0;	## Define desired ORF type 

my $usage = 'perl getorf.pl *.fasta';

while (my $file = shift@ARGV){
	$file =~ s/.fsa//;
	system "echo Working on $file...";
	system "getorf -minsize 300 -find $type -sequence $file.fsa -outseq $file.orf";
}
