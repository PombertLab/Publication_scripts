#!/usr/bin/perl
## Pombert Lab, 2019
my $version = '0.2';
my $name = 'extract_FASTA_from_BLAST_hits.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Extracts sequence from BLAST hits against the NCBI nr database
REQUIREMENTS	NCBI Blast+ (blastdbcmd) - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
		A local copy of the NCBI NR database - ftp://ftp.ncbi.nlm.nih.gov/blast/db/
		
USAGE		extract_FASTA_from_BLAST_hits.pl -db nr -e 1e-05 -b *.blast.6 

OPTIONS:
-db		Database queried
-e		E-value cutoff
-b		BLAST tabular output file(s) (outfmt 6)
OPTIONS
die "$options\n" unless @ARGV;

my @blast;
my $cutoff;
my $db;
GetOptions(
	'b=s@{1,}' => \@blast,
	'e=s' => \$cutoff,
	'db=s' => \$db,
);

while (my $file = shift@blast){
	open IN, "<$file";
	$file =~ s/\..*$//;
	system "echo Working on $file.$cutoff.fasta";
	while (my $line = <IN>){
		chomp $line;
		my @columns = split ("\t", $line);
		my $query = $columns[0];
		my ($hit) = $columns[1] =~ /^gi\|(\d+)/;
		my $evalue = $columns[10];
		if ($evalue <= $cutoff){system "blastdbcmd -entry $hit -db $db -outfmt '%f' >> $file.$cutoff.fasta";}
	}
}
