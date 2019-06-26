#!/usr/bin/perl
## Pombert Lab, 2018
my $name = 'make_datasets.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Creates FASTA datasets from Orthofinder orthogroups output files
EXAMPLE		make_datasets.pl -f mRNAs.fasta -c file.csv

OPTIONS:
-f (--fasta)	FASTA input file (mRNAs)
-c (--csv)	CSV input file (orthogroups from Orthofinder)
OPTIONS
die "$options\n" unless @ARGV;

my $fasta; my $csv;
GetOptions(
	'f|fasta=s' => \$fasta,
	'c|csv=s' => \$csv
);

open RNA, "<$fasta"; open CSV, "<$csv";

## Creating RNA database
my %RNA;
my $locus;
while (my $line = <RNA>){
	chomp $line;
	if ($line =~ /^>(\S+)/){$locus = $1;}
	else{$RNA{$locus} .= $line;}
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
				print OUT ">$sp\n";
				my @nt = unpack("(A60)*", $RNA{$sp});
				while (my $tmp = shift@nt){print OUT "$tmp\n";}
			}
		}
		close OUT;
	}
}