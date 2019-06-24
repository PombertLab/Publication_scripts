#!/usr/bin/perl
## Pombert Lab, 2018
## Create Fasta datasets from Orthofinder Orthgroups output files.

use strict;
use warnings;

die "USAGE = make_datasets.pl mRNAs.fasta file.csv\n" unless @ARGV;

open RNA, "<$ARGV[0]";
open CSV, "<$ARGV[1]";

## Creating RNA database
my %RNA;
my $name;
while (my $line = <RNA>){
	chomp $line;
	if ($line =~ /^>(\S+)/){$name = $1;}
	else{$RNA{$name} .= $line;}
}

## Creating fastas
while (my $line = <CSV>){
	chomp $line;
	if ($line =~ /^\t/){next;}
	else{
		my @OG = split(/, |\t/, $line);
		my $og = shift@OG; print "Working on $og...\n";
		open OUT, ">$og.fasta";
		for my $sp (@OG){
			if ($sp eq ''){next;}
			else{
				(my $tag) = ($sp =~ /^(\w+)_/);
				print OUT ">$tag\n";
				my @nt = unpack("(A60)*", $RNA{$sp});
				while (my $tmp = shift@nt){print OUT "$tmp\n";}
			}
		}
		close OUT;
	}
}