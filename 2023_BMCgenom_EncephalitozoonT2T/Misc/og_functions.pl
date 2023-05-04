#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'og_functions.pl';
my $version = '0.1';
my $updated = '2022-07-08';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Assign functions to OrthoFinder orthogroups based on a reference set

COMMAND		${name} \\
		  -t *.tsv \\
		  -g Orthogroups.tsv \\
		  -o output.txt

-t (--tsv)	Tab-delimited file(s) containing functions to assign to Orthogroups
-g (--og)	Orthogroups.tsv file from OthoFinder
-o (--out)	Output file
-cl		Column containing locus [Default: 2]
-cf		Column containing functions [Default -1]
USAGE
die "\n$usage\n" unless @ARGV;

my @tsv;
my $ortho;
my $output;
my $cl = 2;
my $cf = -1;
GetOptions(
	't|tsv=s@{1,}' => \@tsv,
	'g|og=s' => \$ortho,
	'o|out=s' => \$output,
	'cl=i' => \$cl,
	'cf=i' => $cf
);

## Creating db of functions
my %functions;
while (my $file = shift@tsv){

	open TSV, "<", $file or die "Can't open $file: $!\n";

	while (my $line = <TSV>){
		chomp $line;
		my @data = split("\t", $line);
		my $locus = $data[$cl];
		$locus =~ s/ //g;
		$locus =~ s/\n//g;
		my $function = $data[$cf];
		$functions{$locus} = $function; 
	}

}

## Working on OrthoFinder output file
open OG, "<", $ortho or die "Can't open $ortho: $!\n";
open OUT, ">", $output or die "Can't create $output: $!\n";

while (my $line = <OG>){

	chomp $line;

	unless ($line =~ /^Orthogroup/){

		my @data = split ("\t", $line);
		my $ogroup = shift @data;
		print OUT "### $ogroup\n";

		foreach my $entry (@data){
			if ($entry =~ /\S/){
				$entry =~ s/ //g;
				my @loci = split(",", $entry);
				foreach my $locus (@loci){
					if (exists $functions{$locus}){
						print OUT "$locus\t$functions{$locus}\n";
					}
				}
			}
		}

	}

}