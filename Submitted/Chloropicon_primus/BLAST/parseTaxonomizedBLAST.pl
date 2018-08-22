#!/usr/bin/perl
## Pombert Lab, IIT 2016

use strict;
use warnings;

my $usage = "USAGE = perl parseTaxonomizedBLAST.pl BLAST_output_file(s)";
die "$usage\n" unless @ARGV;

open EUK1, ">Eukaryotes_plants.txt";
open EUK2, ">Eukaryotes_others.txt";
open BAC1, ">Bacteria_cyanobateria.txt";
open BAC2, ">Bacteria_others.txt";
open UNK, ">Unsorted.txt"; ## Keepting track of regex failures

while (my $file = shift@ARGV){
	open IN, "<$file";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(\w+)\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\s+(\S+)\t(\S+)\t(\S+)(.*)\t(Eukaryota|Bacteria)\t(.*)$/){
			my $query = $1;
			my $target = $2;
			my $qstart = $3;
			my $qend = $4;
			my $pident = $5;
			my $length = $6;
			my $bitscore = $7;
			my $evalue = $8;
			my $taxid = $9;
			my $sciname = $10;
			my $kingdom = $11;
			my $blastname = $12;
			print "$blastname\n";
			## Eukaryotes
			if (($kingdom eq 'Eukaryota')&&($blastname eq 'green plants')){print EUK1 "$line\n";}
			elsif (($kingdom eq 'Eukaryota')&&($blastname eq 'green algae')){print EUK1 "$line\n";}
			elsif (($kingdom eq 'Eukaryota')&&($blastname eq 'monocots')){print EUK1 "$line\n";}
			elsif (($kingdom eq 'Eukaryota')&&($blastname eq 'eudicots')){print EUK1 "$line\n";}
			elsif (($kingdom eq 'Eukaryota')&&($blastname eq 'mosses')){print EUK1 "$line\n";}
			elsif (($kingdom eq 'Eukaryota')&&($blastname eq 'club-mosses')){print EUK1 "$line\n";}
			elsif (($kingdom eq 'Eukaryota')&&($blastname eq 'conifers')){print EUK1 "$line\n";}
			elsif (($kingdom eq 'Eukaryota')&&($blastname eq 'flowering plants')){print EUK1 "$line\n";}
			elsif ($kingdom eq 'Eukaryota'){print EUK2 "$line\n";}
			## Bacteria
			elsif (($kingdom eq 'Bacteria')&&($blastname eq 'cyanobacteria')){print BAC1 "$line\n";}
			elsif ($kingdom eq 'Bacteria'){print BAC2 "$line\n";}
			## Keepting track of regex failures
			else {print UNK "$line\n";}
		}
	}
}