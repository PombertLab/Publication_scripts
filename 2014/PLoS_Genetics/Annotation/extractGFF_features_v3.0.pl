#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl extractGFF_features.pl *.gff';

die $usage unless @ARGV;



while (my $gff = shift @ARGV) {

	open IN, "<$gff" or die "cannot open $gff";
	$gff =~ s/\.gff$//;
	open OUT, ">$gff.feat";
	
	print OUT "strand\ttype\tstart\tend\tlocus_tag\n";

	while (my $line = <IN>) {

		chomp $line;
		
		if ($line =~ /^contig-\d+\t\.\tcontig\t\d\t(\d+)\t\.\t\.\t\..*$/) {
			my $contig_length = $1;
			print OUT "length\t$contig_length\n";
		}
		
		elsif ($line =~ /^contig-\d+\taugustus_masked\t(\S+)\t(\d+)\t(\d+)\t\S+\t([+-]).*Name=augustus_masked-contig-(\d+)-abinit-gene-\d+\.(\d+)-mRNA-(\d+);(_AED|Target).*$/) {
			
			my $strandedness = $4;
			my $type = $1;
			my $start = $2;
			my $end = $3;
			my $locus_tag = $5;
			my $locus_tag_num = $6;
			my $locus_tag_num_2 = $7;

			
			print OUT "$strandedness"."\t";
			print OUT "$type"."\t";
			print OUT "$start"."\t";
			print OUT "$end"."\t";
			print OUT "H632_c$locus_tag\.$locus_tag_num"."\n";
		}
	}
	close IN;
}

exit;
