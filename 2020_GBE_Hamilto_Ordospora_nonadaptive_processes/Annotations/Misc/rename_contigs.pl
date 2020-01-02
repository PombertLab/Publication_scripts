#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'rename_contigs.pl';
my $version = '0.2'; 

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Rename/renumber contigs in fasta files
USAGE		rename_contigs.pl -prefix contig -digit 4 -increment 1 -in input.fasta -out output.fasta

OPTIONS:
-prefix		Sequence name prefix [default: contig]
-digit		Number of digits [default: 4]
-increment	Auto-increment value [default: 1]
-in		FASTA input file
-out		Desired output name
OPTIONS
die "$usage\n" unless @ARGV;

my $prefix = 'contig';
my $digit = 4;
my $increment = 1;
my $in;
my $out;
GetOptions(
	'prefix=s' => \$prefix,
	'digit=i' => \$digit,
	'increment=i' => \$increment,
	'in=s' => \$in,
	'out=s' => \$out
);

open IN, "<$in"; open OUT, ">$out";

my $count = 0;
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^>/){
	 	$count += $increment;
		my $num = sprintf("%0${digit}d", $count);
		print OUT ">${prefix}-$num\n";
	}
	else {
		print OUT "$line\n";
	}
}
