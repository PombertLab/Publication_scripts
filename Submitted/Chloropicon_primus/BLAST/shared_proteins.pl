#!/usr/bin/perl
## Pombert Lab IIT, 2016

use strict;
use warnings;

my $usage = 'script.pl min_evalue protein_list *.hits';
die "\nUSAGE = $usage

min_evalue = minimum BLAST Evalue
protein_list = list of proteins queried against databases (one per line)
*.hits = concatenated BLASTP/TBLASTN outputs\n
" unless@ARGV;

my $evalue = shift@ARGV;
my $scifmt = sprintf("%.100f", $evalue);
my $list = shift@ARGV;
my %hash = ();

while (my $file = shift@ARGV){
	open HIT, "<$file";
	open OUT, ">$file.$evalue.shared";
	open LIST, "<$list";
	%hash = ();
	while (my $line = <HIT>){
		chomp $line;
		if ($line =~/^(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/){
			my $locus = $1; my $hit = $2; my $eval = $3;
			my $ev = sprintf("%.100f", $eval);
			if ($ev <= $scifmt){
				$hash{$locus} = $locus;
				print "FOUND: $hash{$locus}\n"; ## Debugging
			}
		}
	}
	while (my $line = <LIST>){
	chomp $line;
		if ($line =~ /^#/){next;}
		elsif ($line =~ /^(\S+)/){
			my $prot = $1;
			if (exists $hash{$prot}){print OUT "$prot\n"}
		}
	}
}