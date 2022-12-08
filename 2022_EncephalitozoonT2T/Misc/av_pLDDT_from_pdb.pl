#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'av_pLDDT_from_pdb.pl';
my $version = '0.1';
my $updated = '2022-10-26';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = << "USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Calculates average pLDDT scores from PDB files from the B-factor field.
		NOTE: Only works for PDB files with pLDDT scores in the B-factor field.

COMMAND		${name} -p *.pdb -o pLDDTs.tsv

OPTIONS:
-p (--pdb)	PDB file(s) to check
-o (--out)	Tab-delimited output file
USAGE
die "\n$usage\n" unless @ARGV;

my @pdbs;
my $outfile;
GetOptions(
	'p|pdb=s@{1,}' => \@pdbs,
	'o|out=s' => \$outfile
);

open OUT, ">", $outfile or die "Can't create $outfile: $!\n";
print OUT "# File\tMolecule\tAverage pLDDT score\n";

while (my $pdb_file = shift@pdbs){

	open PDB, "<", $pdb_file or die "Can't open $pdb_file: $!\n";

	my $sum_pLDDT = 0;
	my $res_counter = 0;
	my $molecule;
	while (my $line = <PDB>){
		chomp $line;
		if ($line =~ /^COMPND.*MOLECULE: (.*);/){
			$molecule = $1;
		}
		if ($line =~ /^ATOM/){
			$res_counter++;
			my $pLDDT = substr($line, -20, 6);
			$pLDDT =~ s/\s//g;
			$sum_pLDDT += $pLDDT;
		}
	}
	my $av_pLDDT = $sum_pLDDT/$res_counter;
	$av_pLDDT = sprintf("%.2f", $av_pLDDT);

	if (!defined $molecule){
		$molecule = 'no entry';
	}

	print OUT $pdb_file."\t".$molecule."\t".$av_pLDDT."\n";

}
