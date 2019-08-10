#!/usr/bin/perl
## Pombert Lab, 2019
my $version = 0.1;
my $name = 'FASTQtoFASTA.pl';

use strict; use warnings;

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Converts FASTQ files to FASTA format (without quality scores)
USAGE		FASTQtoFASTA.pl *.fastq
OPTIONS
die "$usage\n" unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.fastq$//; $file =~ s/fq$//;
	open OUT, ">$file.fasta";
	my $count = 0;
	while (my $line = <IN>){
		chomp $line;
		if ($count == 0){
			$line =~ s/^@/>/;
			print OUT "$line\n";
			$count++;
		}
		elsif ($count == 1){
			print OUT "$line\n";
			$count++;
		}
		elsif ($count == 2){$count++;}
		elsif ($count == 3){$count = 0;}
	}
}