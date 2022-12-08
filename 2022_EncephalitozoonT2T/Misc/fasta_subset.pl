#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'fasta_subset.pl';
my $version = '0.1';
my $updated = '2022-07-07';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = << "USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Creates FASTA subsets from a list in TSV format 

COMMAND		${name} \\
		  -f file.fasta \\
		  -t file.tsv \\
		  -o subset.fasta

-f (--fasta)	FASTA file(s) containing the sequences to extract
-t (--tsv)	TSV file containing list of desired sequences
-c (--column)	Column containing the desired sequences [Default: 3]
-o (--output)	Output file containing desired subset of sequences
USAGE
die "\n$usage\n" unless @ARGV;

my @fasta;
my $tsv;
my $column = 3;
my $output;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	't|tsv=s' => \$tsv,
	'c|column=i' => \$column,
	'o|output=s' => \$output
);

### Creating a database of sequences from the FASTA file(s)
my %sequences;
while (my $fasta = shift @fasta){

	open FASTA, "<", $fasta or die "Can't open $fasta: $!\n";
	my $sequence;

	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			$sequence = $1;
			# print "$sequence\n";
		}
		else {
			$sequences{$sequence} .= $line;
		}
	}
}

# exit;

### Parsing the TSV file and extracting the sequences from the db
open TSV, "<", $tsv or die "Can't open $tsv: $!\n";
open OUT, ">", $output or die "Can't create $output: $!\n";

while (my $line = <TSV>){

	chomp $line;
	unless ($line =~ /^#/){
		my @data = split ("\t", $line);
		my $seq = $data[$column - 1];
		if (exists $sequences{$seq}){
			print OUT ">$seq\n";
			my @tmp = unpack ("(A60)*", $sequences{$seq});
			while (my $tmp = shift@tmp){
				print OUT "$tmp\n";
			}
		}
		else {
			print STDERR "[Error]: Sequence $seq not found in FASTA file(s) provided...\n";
		}
	}

}