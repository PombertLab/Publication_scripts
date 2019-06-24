#!/usr/bin/perl

use strict;
use warnings;

die "\nadd_prefixes.pl NEWPREFIX *.embl\n\n" unless @ARGV;

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