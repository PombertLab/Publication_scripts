#!/usr/bin/perl
## Pombert Lab, 2019
my $name = 'low_gc.pl';
my $version = 0.2;

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Returns genomic regions with GC content lower/higher than specifed value
USAGE		low_gc.pl -max 50 -gc genome.gc

OPTIONS:
-max		GC cutoff <= value [Default: 50]
-gc		GC content file (generated with GC_content_to_Circos.pl) [Default: genome.gc]
OPTIONS
die "$options\n" unless @ARGV;

my $gc = 50;
my $file = 'genome.gc';
GetOptions(
	'max=i' => \$gc,
	'gc=s' => \$file
);


open IN, "<$file";
open LOW, ">$file.$gc.low_gc"; ## <= specified GC cutoff
open HIGH, ">$file.$gc.high_gc"; ## > specified GC cutoff
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^#/){next;}
	else{
		my @cols = split(" ", $line);
		if ($cols[3] <= $gc){print LOW "$line\n";}
		else {print HIGH "$line\n";}
	}
}
