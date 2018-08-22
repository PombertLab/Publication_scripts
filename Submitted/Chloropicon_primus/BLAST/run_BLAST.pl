#!/usr/bin/perl
## Pombert Lab, IIT 2017

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## User 
my $dbdir = "/media/Data_1/jpombert/CCMP1205/BLAST/Table_S1/DB/"; ## Path to databases

## Usage definition
my $usage = "\nUSAGE = run_BLAST.pl [options]\n
EXAMPLE: run_BLAST.pl -type blastp -e 1e-15 -db databases.txt -fa *.fasta -t 10\n";
my $hint = "Type run_TBLASTN.pl -h (--help) for list of options\n";
die "$usage\n$hint\n" unless@ARGV;

## Defining options
my $options = <<'END_OPTIONS';
OPTIONS:
-h (--help)	Displays this list of options
-type		Type of blast; blastp, tblastn, blastx, blastn [default: blastp]
-t (--threads)	Number of CPU threads to use [default: 10]
-e (--evalue)	Minimum evalue [default: 1e-05]
-fa (--fasta)	Fasta files to investigate
-db		List of databases to query against, one per line
END_OPTIONS

my $help = '';
my $type = 'blastp';
my $threads = '10';
my $evalue = 1e-05;
my @fasta ;
my $dblist = '';

GetOptions(
	'h|help' => \$help,
	'type=s' => \$type,
	't|threads=i' => \$threads,
	'e|evalue=s' => \$evalue,
	'fa|fasta=s@{1,}' => \@fasta,
	'db=s' => \$dblist
);

if ($help){die "$usage\n$options";}

while (my $file = shift@fasta){
	open IN, "<$dblist";
	while (my $db= <IN>){
		chomp $db;
		system "echo Performing $type using input file $file against database $db..."; 
		system "$type -num_threads $threads -query $file -db $dbdir"."$db -evalue $evalue -culling_limit 1 -outfmt 6 -out $file.vs.$db.$evalue.$type";
	}
}

