#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl RNAmmer_to_EMBL.pl *.gff3';
die "$usage\n" unless @ARGV;

my %RNA = (
'28s_rRNA' => 'large subunit ribosomal RNA',
'23s_rRNA' => 'large subunit ribosomal RNA',
'18s_rRNA' => 'small subunit ribosomal RNA',
'16s_rRNA' => 'small subunit ribosomal RNA',
'8s_rRNA' => '5S ribosomal RNA',
'5s_rRNA' => '5S ribosomal RNA',
);

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.gff//;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(\w+)\s+RNAMMER\s+rRNA\t(\d+)\t(\d+)\t\S+\t(\+|\-).*Note=(\w+)/){
			my $chromosome = $1;
			my $start = $2;
			my $end = $3;
			my $strand = $4;
			my $product = $5;
			open OUT, ">>$chromosome.RNAmmer.embl";
			if ($strand eq '+'){
				print OUT "FT   gene             $start..$end\n";
				print OUT "FT   rRNA             $start..$end\n";
				print OUT "FT                   /product=\"$RNA{$product}\"\n";
			}
			elsif ($strand eq '-'){
				print OUT "FT   gene             complement($end..$start)\n";
				print OUT "FT   rRNA             complement($end..$start)\n";
				print OUT "FT                   /product=\"$RNA{$product}\"\n";
			}
		}
	}
}
