#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl rename_proteins.pl *.fasta';

die $usage unless @ARGV;



while (my $fasta = shift @ARGV) {

	open IN, "<$fasta" or die "cannot open $fasta";
	$fasta =~ s/.maker.non_overlapping_ab_initio.proteins.fasta$//;
	open OUT, ">$fasta.proteins.fsa";

	while (my $line = <IN>) {

		chomp $line;

		if ($line =~ />augustus_masked-contig-(\d+)-abinit-gene-(\d+).(\d+)-mRNA-(\d+).*$/) {
			print OUT ">c$1\-$2\.$3\-$4\n";
		}
		
		else {
			print OUT "$line\n";
		}
	}
	close IN;
}

exit;
