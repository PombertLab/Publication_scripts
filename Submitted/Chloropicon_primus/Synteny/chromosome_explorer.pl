#!/usr/bin/perl
## Pombert Lab, 2017
## v0.1

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = "USAGE = chromosome_explorer.pl -q CCMP.list -b Chlorella.blastp.6 -s Chlorella.list -o CCMPvsChlorella.chromosomes";
die "\n$usage\n\n" unless @ARGV;

my $options = <<"OPTIONS";

USAGE = chromosome_explorer.pl -q CCMP.list -b Chlorella.blastp.6 -s Chlorella.list -o CCMPvsChlorella.chromosomes

-h (--help)		Display this help message
-q (--query)		Query input list ## Must correspond to query (qseqid) in blast output
-b (--blast)		BLASTP output in format 6
-s (--subject)		Subject input list ## Must correspond to subject (sseqid) in blast output
-o (--output)		Output file ## Writes chromosome/contigs metrics

OPTIONS

my $help;
my $query;
my $blast;
my $subject;
my $output;

GetOptions(
	'h|help' => \$help,
	'q|query=s' => \$query,
	'b|blast=s' => \$blast,
	's|subject=s' => \$subject,
	'o|output=s' => \$output
);
die "$options" if $help;

## Parsing BLAST
open BLAST, "<$blast" or die "Can't open blast file : $blast\n";
my %hits;
while (my $line = <BLAST>){
	chomp $line;
	if ($line =~ /^(\S+)\t(\S+)/){
		my $query = $1;
		my $hit = $2;
		if (!exists $hits{$query}){$hits{$query} = $hit;}
	}
}

## Creating target genome DB
open LIST, "<$subject" or die "Can't open subject file : $subject\n";
my %ortho;
while (my $line = <LIST>){
	chomp $line;
	if ($line =~ /^(\S+)\t(\S+)\s+(\d+)\s+(\d+)\t([+-])\t(\d+)/){
		my $locus = $1; my $contig = $2;
		$ortho{$locus} = $contig;
	}
}

## Working on our CCMP.genome.list
open GENOME, "<$query" or die "Can't open query file : $query\n";
open OUT, ">>$output" or die "Can't open $output\n";
open OUT2, ">>$output.metrics";
my %hash;
my $prot = 0; my $hit = 0;
while (my $line = <GENOME>){
	chomp $line;
	if ($line =~ /^(\S+)\t(\S+)\s+(\d+)\s+(\d+)\t([+-])\t(\d+)/){
		my $query = $1; my $chromosome = $2; $prot++;
		if (!exists $hits{$query}){next;}
		else{
			$hit++;
			if (!exists $hash{$ortho{$hits{$query}}}){$hash{$ortho{$hits{$query}}} = 1;}
			else{$hash{$ortho{$hits{$query}}} += 1;}
		}
	}
}
print OUT "## $query\n";
foreach my $name (sort { $hash{$b} <=> $hash{$a} } keys %hash) {
	printf OUT "%-8s %s\n", $name, $hash{$name};
}
print OUT2 "$query\tproteins = $prot\t";
print OUT2 "blast hits = $hit\n";
close OUT;
close OUT2;