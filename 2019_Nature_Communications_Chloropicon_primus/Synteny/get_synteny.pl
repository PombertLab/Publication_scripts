#!/usr/bin/perl
## Pombert Lab, 2017
## v0.2

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = "USAGE = get_synteny.pl -q CCMP.list -b Chlorella.blastp.6 -s Chlorella.list -g 10 -o CCMPvsChlorella";
die "\n$usage\n\n" unless @ARGV;

my $options = <<"OPTIONS";

USAGE = get_synteny.pl -q CCMP.list -b Chlorella.blastp.6 -l Chlorella.list -gap 10

-h (--help)		Display this help message
-q (--query)		Query input list ## Must correspond to query (qseqid) in blast output
-b (--blast)		BLASTP output in format 6
-s (--subject)		Subject input list ## Must correspond to subject (sseqid) in blast output
-g (--gap)		Space allowed between adjacent genes [Default: 0]
-o (--output)		Output file ## Writes detected gene pairs

OPTIONS

my $help;
my $query;
my $blast;
my $subject;
my $gap = 50;
my $output;

GetOptions(
	'h|help' => \$help,
	'q|query=s' => \$query,
	'b|blast=s' => \$blast,
	's|subject=s' => \$subject,
	'g|gap=i' => \$gap,
	'o|output=s' => \$output
);
die "$options" if $help;

## Parsing BLAST
open BLAST, "<$blast" or die "Can't open blast file : $blast\n";
my %hits;
while (my $line = <BLAST>){
	chomp $line;
	if ($line =~ /^(\S+)\t(\S+)/){
		my $query = $1;
		my $hit = $2;
		if (!exists $hits{$query}){$hits{$query} = $hit;}
	}
}

## Creating target genome DB
open LIST, "<$subject" or die "Can't open subject file : $subject\n";
my %ortho;
while (my $line = <LIST>){
	chomp $line;
	if ($line =~ /^(\S+)\t(\S+)\s+(\d+)\s+(\d+)\t([+-])\t(\d+)/){
		my $locus = $1; my $contig = $2;	my $start = $3; my $end = $4; my $strand = $5; my $pos = $6;
		$ortho{$locus}[0] = $contig;	#print "contig = $ortho{$locus}[0]\n"; ## Contig of chromosome
		$ortho{$locus}[1] = $strand;	#print "strand = $ortho{$locus}[1]\n"; ## Strandedness
		$ortho{$locus}[2] = $pos;	#print "pos = $ortho{$locus}[2]\n";	## Position on the contig/chromosome 
	}
}

## Working on our CCMP.genome.list
open GENOME, "<$query" or die "Can't open query file : $subject\n";
open OUT, ">$output.gap$gap.clusters";
print OUT "##Query_QA\tStrand_QA\tQuery_QB\tStrand_QB\tQuery_contig\tSubject_SA\tStrand_SA\tSubject_SB\tStrand_SB\tSubject_contig\tRelative_orientation\tGenes_in-between\n"; 
my $previous_hit;
my $previous_strand;
while (my $line = <GENOME>){
	chomp $line;
	if ($line =~ /^(\S+)\t(\S+)\s+(\d+)\s+(\d+)\t([+-])\t(\d+)/){
		my $query = $1; my $chromosome = $2; my $start = $3; my $end = $4; my $strand = $5; my $pos =  $6;
		if (!exists $hits{$query}){next;}
		elsif ( (($strand eq '-') && ($ortho{$hits{$query}}[1] eq '-')) || (($strand eq '+') && ($ortho{$hits{$query}}[1] eq '+')) ){	## Testing for same direction
			if (!defined $previous_hit){$previous_hit = $query; $previous_strand = $strand;}
			else{ ## Must check for same chromosome [0], relative strandedness [1] and positions [2]
				if (($ortho{$hits{$query}}[0] eq $ortho{$hits{$previous_hit}}[0])  && ($previous_strand eq $ortho{$hits{$previous_hit}}[1]) && (($ortho{$hits{$query}}[2] > $ortho{$hits{$previous_hit}}[2]) && ($ortho{$hits{$query}}[2] <= ($ortho{$hits{$previous_hit}}[2] + 1 + $gap)))){
					my $sep = ($ortho{$hits{$query}}[2] - $ortho{$hits{$previous_hit}}[2]  - 1);
					print OUT "$previous_hit\t$previous_strand\t$query\t$strand\t$chromosome\t$hits{$previous_hit}\t$ortho{$hits{$previous_hit}}[1]\t$hits{$query}\t$ortho{$hits{$query}}[1]\t$ortho{$hits{$query}}[0]\tFORWARD\t$sep\n";
				}
				$previous_hit = $query; $previous_strand = $strand;
			}
		}
		elsif ( (($strand eq '+') && ($ortho{$hits{$query}}[1] eq '-')) || (($strand eq '-') && ($ortho{$hits{$query}}[1] eq '+')) ){	## Testing for opposite directions
			if (!defined $previous_hit){$previous_hit = $query; $previous_strand = $strand;}
			else{ ## Must check for same chromosome [0], relative strandedness [1] and positions [2]
				if (($ortho{$hits{$query}}[0] eq $ortho{$hits{$previous_hit}}[0])  && ($previous_strand ne $ortho{$hits{$previous_hit}}[1]) && (($ortho{$hits{$query}}[2] < $ortho{$hits{$previous_hit}}[2]) && ($ortho{$hits{$query}}[2] >= ($ortho{$hits{$previous_hit}}[2] - 1 - $gap)))){
					my $sep = ($ortho{$hits{$previous_hit}}[2]  - $ortho{$hits{$query}}[2]  - 1);
					print OUT "$previous_hit\t$previous_strand\t$query\t$strand\t$chromosome\t$hits{$previous_hit}\t$ortho{$hits{$previous_hit}}[1]\t$hits{$query}\t$ortho{$hits{$query}}[1]\t$ortho{$hits{$query}}[0]\tREVERSE\t$sep\n";
				}
				$previous_hit = $query; $previous_strand = $strand;
			}
		}
	}
}
