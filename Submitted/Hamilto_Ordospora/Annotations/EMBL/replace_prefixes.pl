#!/usr/bin/perl
## Pombert Lab, IIT, 2018
my $name = 'replace_prefixes.pl';
my $version = '0.1';

use strict; use warnings;

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Replaces locus_tag prefixes in EMBL files; useful if temporary prefixes were affixed while waiting for the NCBI ones...
USAGE		replace_prefixes.pl NEWPREFIX *.embl
OPTIONS
die "$usage\n" unless @ARGV;

my $prefix = shift@ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.tmp";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(FT\s+\/locus_tag=")(\w+)(_\S+)$/){
			print OUT "$1$prefix$3\n";
		}
		else{print OUT "$line\n";}
	}
	close IN; close OUT;
	system "mv $file.tmp $file";
}