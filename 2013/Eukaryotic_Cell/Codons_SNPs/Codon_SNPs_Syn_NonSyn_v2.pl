#!/usr/bin/perl

## Query each SNP-containing codon against the hash to verify if the resulting amino acid is the same.
## Write synonymous changes to .syn, non-synonymous changes to non-syn.
## The codon hash table can be modified to reflect modifications to the genetic code. 

use strict;
use warnings;

my $usage = 'perl Codon_SNPs_Syn_NonSyn.pl all_codon_snps.snps';

die $usage unless @ARGV;

my %aa = (
'tca'=>'S','tcc'=>'S','tcg'=>'S','tct'=>'S',
'ttc'=>'F','ttt'=>'F',
'tta'=>'L','ttg'=>'L',
'tac'=>'Y','tat'=>'Y',
'taa'=>'_','tag'=>'_','tga'=>'_',
'tgc'=>'C','tgt'=>'C',
'tgg'=>'W',
'cta'=>'L','ctc'=>'L','ctg'=>'L','ctt'=>'L',
'cca'=>'P','ccc'=>'P','ccg'=>'P','cct'=>'P',
'cac'=>'H','cat'=>'H',
'caa'=>'Q','cag'=>'Q',
'cga'=>'R','cgc'=>'R','cgg'=>'R','cgt'=>'R',
'ata'=>'I','atc'=>'I','att'=>'I',
'atg'=>'M',
'aca'=>'T','acc'=>'T','acg'=>'T','act'=>'T',
'aac'=>'N','aat'=>'N',
'aaa'=>'K','aag'=>'K',
'agc'=>'S','agt'=>'S',
'aga'=>'R','agg'=>'R',
'gta'=>'V','gtc'=>'V','gtg'=>'V','gtt'=>'V',
'gca'=>'A','gcc'=>'A','gcg'=>'A','gct'=>'A',
'gac'=>'D','gat'=>'D',
'gaa'=>'E','gag'=>'E',
'gga'=>'G','ggc'=>'G','ggg'=>'G','ggt'=>'G');


while (my $snps = shift @ARGV) {

	open IN, "<$snps" or die "cannot open $snps";
	$snps =~ s/.snps//;
	open SYN, ">$snps.syn";
	open NONSYN, ">$snps.nonsyn";
	
		while (my $line = <IN>) {

			chomp $line;
		
			if ($line =~ /^(ECU\d+_\S+)\t(\d+)\t(\d+)\t(\w{3})\t(\w{3})\t(\w{3})/) {
			
			my $gene = $1;
			my $position = $2;
			my $length =	$3;
			my $strain1 = $4;
			my $strain2 = $5;
			my $strain3 = $6;
			
				if (($aa{$strain1} eq $aa{$strain2}) && ($aa{$strain1} eq $aa{$strain3})) {
				
				print SYN "$gene\t$position\t$length\t$aa{$strain1}\t$aa{$strain2}\t$aa{$strain3}\n";
				
				}
				
				else {
				print NONSYN 	"$gene\t$position\t$length\t$aa{$strain1}\t$aa{$strain2}\t$aa{$strain3}\n";
				}
			
			}	
		}

}

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)