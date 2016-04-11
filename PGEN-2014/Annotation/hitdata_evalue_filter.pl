#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl hitdata_evalue_filter.pl hitdata.txt';

die $usage unless @ARGV;

while (my $txt = shift @ARGV) {

	open IN, "<$txt" or die "cannot open $txt";
	open OUT, ">$txt.evalues";
	
		while (my $line = <IN>) {

			chomp $line;
			
			if ($line =~ /^Q#\d+ - >(H632_c\d+.\d+)\t\S+\t\d{1,}\t\d{1,}\t\d{1,}\t0\.\d{1,}/) {
				next
			}
			elsif ($line =~ /^Q#\d+ - >(H632_c\d+.\d+)\t\S+\t\d{1,}\t\d{1,}\t\d{1,}\t(0)/) {
				print OUT "$line\n";
			}
			elsif ($line =~ /^Q#\d+ - >(H632_c\d+.\d+)\t\S+\t\d{1,}\t\d{1,}\t\d{1,}\t\d{1,}.\d{1,}e-(\d{1,})/){
				my $evalue = $2;
				
				if ($evalue >= 30){
				print OUT "$line\n";
				}
			}
		}
}


