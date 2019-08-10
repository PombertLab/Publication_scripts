#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl whitespaces_to_tab.pl file.txt';

die $usage unless @ARGV;

while (my $txt = shift @ARGV) {

	open IN, "<$txt" or die "cannot open $txt";
	open OUT, ">$txt.tabs";
	
		while (my $line = <IN>) {

			chomp $line;
			
			if ($line =~ /^H632/) {
				
				$line =~ s/\s{1,}/\t/g;
				print OUT "$line\n";
			}
		}
}


