#!/usr/bin/perl
## Pombert Lab, 2018
my $name = 'Blast_to_Circos.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generates highlights from BLAST hits that can be used as tracks for Circos
EXAMPLE		Blast_to_Circos.pl -l genome.list -b file.blastp.6 -c red -e 1e-10

OPTIONS:
-l	List of locus_tags and their positions for Circos ## See genome.list
-b	BLAST output (outfmt 6 format) of sequences from genome.list as queries
-c	Desired color for the Circos highlights
-e	Minimum BLAST E-value to plot [Default: 1e-10]
OPTIONS
die "$options\n" unless @ARGV;

my $list;
my @blast;
my $evalue = '1e-10';
my $color = 'red';
GetOptions(
	'l=s' => \$list,
	'b=s@{1,}' => \@blast,
	'e=s' => \$evalue,
	'c=s' => \$color
);

## Parsing blast hits
my %hits;
my @species;
while (my $blast = shift@blast){
	open IN, "<$blast";
	$blast =~ s/.blastp.6//;
	push (@species, $blast);
	while (my $line = <IN>){
		chomp $line;
		my @columns = split ("\t", $line); ## query = $columns[0], evalue = $columns[10]
		if (exists $hits{$blast}{$columns[0]}){next;}
		elsif ($columns[10] <= $evalue){
			$hits{$blast}{$columns[0]} = $columns[10];
		}
	}
}
## Working on list(s)
open LIST, "<$list";
while (my $sp = shift@species){
	open OUT, ">$sp.$evalue.list";
	my %hash;
	my $key;
	my $pnum;
	my @keys;
	while (my $line = <LIST>){
		chomp $line;
		if ($line =~ /^(A3770_\d+p)(\d+)\t(\S+)\t(\d+)\t(\d+)/){
			my $prefix = $1; my $num =$2; my $contig = $3;
			my $start = $4; my $end = $5; my $name = $prefix.$num;
			if (exists $hits{$sp}{$name}){
				if (!defined $pnum){
					$pnum = $num;
					$key = $start;
					$hash{$contig}{$key} = $end;
					push (@keys, $contig.$key);
				}
				elsif ($num == $pnum + 10){
					$hash{$contig}{$key} = $end;
					$pnum = $num;
				}
				else{
					$key = $start; $pnum = $num;
					$hash{$contig}{$key} = $end;
					push (@keys, $contig.$key);
				}
			}
		}
	}
	while (my $tile = shift@keys){
		if ($tile =~ /(Chromosome_\d{2})(\d+)/){
			print OUT "$1 $2 $hash{$1}{$2} fill_color=$color,stroke_color=$color\n";
		}
	}
	close OUT;
}
