#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl parse_blast.pl blastoutput_by_species'; ## Use BLAST tabular output format (-outfmt 6)

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.parsed";
	print OUT "Trinity_RNA\tProtein_ID\tE-value\n";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(a\d+\;\d+)\tjgi\|(\w+)\|(\d+)\|\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t(\S+)/){
			my $rna = $1;
			my $species = $2;
			my $proteinid = $3;
			my $evalue = $4;
			
			print OUT "$rna\t$proteinid\t$evalue\n";
		
		}
	}
}



