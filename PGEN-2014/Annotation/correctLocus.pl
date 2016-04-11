#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.test//;
	$file =~ s/.tbl4//;
	open OUT, ">$file.tbl";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(.*H632_c\d+)\.(\d+.*)$/){
			print OUT "$1p$2\n";
		}else{
			print OUT "$line\n";
		}
	}
}