#!/usr/bin/perl
## Pombert Lab, IIT, 2018
my $name = 'add_taxinfo.pl';
my $version = '0.1';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Adds taxonomic info to FASTA headers; useful for TBL2ASN
USAGE		add_taxinfo.pl -sp 'Hamiltosporidium magnivora' -is BE-OM-2 -fa *.fsa
NOTE		Quick script, most info is hard-coded in the print OUT lines...

OPTIONS:
-sp	## species
-is	## isolate
-fa	## FASTA files
OPTIONS
die "$usage\n" unless @ARGV;

my $species;
my $isolate;
my @fasta;
GetOptions(
	'sp=s' => \$species,
	'is=s' => \$isolate,
	'fa=s@{1,}' => \@fasta
);

while (my $file = shift@fasta){
	open IN, "<$file";
	open OUT, ">$file.tmp";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>/){
			print OUT "$line ";
			print OUT '[organism='."$species".']';
			print OUT '[isolate='."$isolate".']';
			print OUT '[lineage=cellular organisms; Eukaryota; Fungi/Metazoa group; Fungi; Microsporidia;][gcode=1]';
			print OUT '[tech=WGS]';
			print OUT '[moltype=DNA]';
			print OUT '[topology=linear]';
			print OUT "\n";
		}
		else{
			print OUT "$line\n";
		}
	}
	close IN;
	close OUT;
	system "mv $file.tmp $file";
}