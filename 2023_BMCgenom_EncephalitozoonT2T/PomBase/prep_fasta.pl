#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'prep_fasta.pl';
my $version = '0.1';
my $updated = '2022-07-29';

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);

my $usage =<< "USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Concatenates FASTA files from PomBase (adds newlines and removes *)

COMMAND		${name} -f *.fasta -out pombase.fasta

-f (--fasta)	FASTA input file from PomBase
-o (--out)	Concatenate FASTA output file
USAGE
die "\n$usage\n" unless @ARGV;

my @fastas;
my $outfile;
GetOptions(
	'f|fasta=s@{1,}' => \@fastas,
	'o|out=s' => \$outfile
);

open OUT, ">", "$outfile" or die "Can't create $outfile: $!\n";

while (my $fasta = shift@fastas){

	my $basename = fileparse($fasta);
	my ($prefix) = $basename =~ /^(.*)\.\w+$/; 

	open FASTA, "<", $fasta or die "Can't open $fasta: $!\n";

	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>/){
			print OUT ">$prefix\n";
		}
		else {
			$line =~ s/\*//g;
			print OUT "$line\n";
		}
	}

}