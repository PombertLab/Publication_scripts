#!/usr/bin/perl
## Pombert Lab, 2017
## blastn -query BEOM2.contigs.fasta -db DB/BEOM2 -outfmt '6 qseqid sseqid qlen slen evalue nident pident' -evalue 1e-50 -out BEOM2.blastn.6
my $name = 'parse_BLAST_selftest.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Parses the output of BLASTN searches of a genome assembly against itself to identify possible duplications
USAGE		parse_BLAST_selftest.pl -p 80 -b BEOM2.blastn.6

OPTIONS:
-p	## Minimum percentage identity [default: 80]
-b	## BLASTN input file; must be performed with the -outfmt '6 qseqid sseqid qlen slen evalue nident pident' option 
OPTIONS
die "$usage\n" unless @ARGV;

my $percent = 80;
my $blastn;
GetOptions(
	'p=i' => \$percent,
	'b=s' => \$blastn
);
my $min = $percent;

## Working on BLASTN file
open IN, "<$blastn";
open OUT, ">$blastn.$min.duplications";
print OUT "Query\tSubject\t# of identities\t% of identities\t% of query\t% of target\tRedundant bp approx\n";

my %contigs;
my $redundancy = 0;
while (my $line = <IN>){
	chomp $line;
	my @columns = split("\t", $line);
	## qseqid0 sseqid1 qlen2 slen3 evalue4 nident5 pident6
	if ($columns[0] eq $columns[1]){next;} ## Ignoring matches between contigs and themselves
	if (exists $contigs{$columns[1]}{$columns[0]}){next;} ## Disregarding query - subject pairs seen previously (removing duplicates)
	else{
		$contigs{$columns[0]}{$columns[1]} = 'seen'; ## Keeping track of query - subject pairs
		my $qpc = sprintf("%.2f", ($columns[5]/$columns[2]*100)); ## Calculating percentage of nucleotides in query that are identical with the subject
		my $spc = sprintf("%.2f", ($columns[5]/$columns[3]*100)); ## Calculating percentage of nucleotides in subject that are identical with the query
		if (($qpc >= $min) || ($spc >= $min)){
			if (($columns[6] >= $min) && ($qpc >= $min)){print OUT "$columns[0]\t$columns[1]\t$columns[5]\t$columns[6]\t$qpc\t$spc\t$columns[2]\n"; $redundancy += $columns[2];}
			elsif (($columns[6] >= $min) && ($spc >= $min)){print OUT "$columns[0]\t$columns[1]\t$columns[5]\t$columns[6]\t$qpc\t$spc\t$columns[2]\n"; $redundancy += $columns[3];}
		}
	}
}
open LOG, ">>Redundancy.txt";
print LOG "Approx. $redundancy bp in $ARGV[1] in are redundant at $min %\n";
