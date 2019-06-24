#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl splitMakerGFF.pl *.gff';
die "$usage\n" unless @ARGV;

while(my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.gff$//;
	open OUT1, ">$file.augustus.gff";
	open OUT2, ">$file.genemark.gff";
	open OUT3, ">$file.repeats.gff";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^\S+\taugustus_masked/){
			print OUT1 "$line\n";
		}
		elsif ($line =~ /^\S+\tgenemark/){
			print OUT2 "$line\n";
		}
		elsif ($line =~ /^\S+\trepeatmasker/){
			print OUT3 "$line\n";
		}
	}
}