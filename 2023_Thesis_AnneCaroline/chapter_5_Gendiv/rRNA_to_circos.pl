#!/usr/bin/env perl
## Pombert Lab, Illinois Tech, 2023

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $name = 'rRNA_to_circos.pl';
my $version = '0.1';
my $updated = '2023-04-06';

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Creates a Circos highlight output file

NOTE		Color names: -- http://circos.ca/tutorials/lessons/configuration/files/images

COMMAND		${name} -g gile.gbff -o output.circos -c blue

OPTIONS:
-g (--genbank)	GenBank file to parse (in gb/gbk/gbff format)
-o (--output)	Desired output file
-c (--color)	Desired color for highlights [Default: blue]
OPTIONS

die "\n$usage\n" unless @ARGV;

my $gbff;
my $circos;
my $color = 'blue';
GetOptions(
	'g|genbank=s' => \$gbff,
	'o|output=s' => \$circos,
	'c|color=s' => \$color
);


open GBK, "<", $gbff or die $!;
open CIRCOS, ">", $circos or die $!;

print CIRCOS "#chr START END COLOR\n";

my $locus;
while (my $line = <GBK>){

	## Grabbing locus from version number
	if ($line =~ /^VERSION\s+(\S+)/){
		$locus = $1;
	}

	## Checking for rRNAs on - strand
	elsif ($line =~ /^\s+rRNA.*\((\d+)..(\d+)\)/){
		my $start = $1;
		my $end = $2;
		print CIRCOS "$locus $start $end fill_color=$color\n";
	}

	## Checking for rRNAs on + strand
	elsif ($line =~ /^\s+rRNA\s+(\d+)..(\d+)/){
		my $start = $1;
		my $end = $2;
		print CIRCOS "$locus $start $end fill_color=$color\n";
	}
}