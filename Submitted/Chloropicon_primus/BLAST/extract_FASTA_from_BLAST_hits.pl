#!/usr/bin/perl

use strict;
use warnings;

my $usage = "USAGE = perl extract_FASTA_from_BLAST_hits.pl db (e.g. nr) cutoff (e.g. 1e-05) blast_hit_files";
die "$usage\n" unless@ARGV;

my $db = shift@ARGV;
my $cutoff = shift@ARGV;
my $min = printf("%.100f", $cutoff); ## Reformatting evalues in decimal numbers

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/\..*$//;
	system "echo Working on $file.$cutoff.fasta";
	system "cat $file.fsa > $file.$cutoff.fasta"; ## Adding the fasta from the query to the top of the file
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(\w+)\tgi\|(\d+)\S+\t\d+\t\d+\t\S+\t\d+\s+\S+\s+(\S+)/){
			my $query = $1;
			my $hit = $2;
			my $evalue = $3;
			my $value = printf("%.100f", $evalue);
			if ($value <= $min){
				system "blastdbcmd -entry $hit -db nr -outfmt '%f' >> $file.$cutoff.fasta";
			}
		}
	}
}