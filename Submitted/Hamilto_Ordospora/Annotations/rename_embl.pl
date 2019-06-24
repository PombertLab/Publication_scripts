#!/usr/bin/perl

use strict;
use warnings;

my $usage = "USAGE = rename.embl.pl PREFIX *.embl";
die "\n$usage\n\n" unless @ARGV;

my $prefix = shift@ARGV;
my @files = @ARGV;
while (my $embl = shift@ARGV){
	open IN, "<$embl";
	open OUT, ">$embl.2";
	my $locus_tag = 0;
	my $contig = $embl; $contig =~ s/.embl$//; $contig =~ s/contig-//i;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^FT\s+gene/){
			$locus_tag += 10; $locus_tag  = sprintf("%04d", $locus_tag); 
			print OUT "$line\n";
		}
		elsif ($line =~ /FT\s+\/locus_tag=/){print OUT "FT                   /locus_tag=\"${prefix}_${contig}p$locus_tag\"\n";}
		else{print OUT "$line\n";}
	}
}
system "rm *.embl";
while (my $file = shift@files){system "mv $file.2 $file";}