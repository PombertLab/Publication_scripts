#!/usr/bin/perl

use strict;
use warnings;

my $usage = "USAGE = filterGC.pl minGC maxGC *.fasta";
die "\n$usage\n\n" unless @ARGV;

my $min = shift@ARGV;
my $max = shift@ARGV;

while (my $fasta = shift@ARGV){
	open IN, "<$fasta";
	open MIN, ">$fasta.mingc";
	open MED, ">$fasta.medgc";
	open MAX, ">$fasta.maxgc";
	my %sequences;
	my $contig;
	my @contigs;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(.*)$/){
			push (@contigs, $1);
			$contig = $1;
		}
		else{
			$sequences{$contig} .= $line;
		}
	}
	while (my $cg = shift@contigs){
		my $gc = 0;
		my $lenght = length($sequences{$cg});
		for (my $i = 0; $i <= $lenght - 1;  $i++){
			my $nt = substr($sequences{$cg}, $i, 1);
			if ($nt eq 'G' || $nt eq 'g' || $nt eq 'C' || $nt eq 'c'){$gc++;}
		}
		my $percent = ($gc/$lenght)*100;
		if ($percent < $min){
			print "UNDER\t$cg\tlenght = $lenght\tGC percent = $percent\n";
			print MIN ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MIN "$sq\n";}
		}
		elsif ($percent > $max){
			print "OVER\t$cg\tlenght = $lenght\tGC percent = $percent\n";
			print MAX ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MAX "$sq\n";}
		}
		else{
			print "BETWEEN\t$cg\tlenght = $lenght\tGC percent = $percent\n";
			print MED ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MED "$sq\n";}
		}
	}
}