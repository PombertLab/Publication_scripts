#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = "EXAMPLE = parse_tcdb.pl -t tcdb.fasta -e 1e-10 -b *.blastp.6 -top 5

-t	TCDB fasta file
-e	Evalue cutoff
-b	BLAST ouput [in outfmt 6 format]
-top	Keep top X hits only [default = 1]
";
die "\n$usage\n" unless @ARGV;

my $tcdb;
my $evalue = '1e-10';
my @blast;
my $top = 1;

GetOptions(
	't=s' => \$tcdb,
	'e=s' => \$evalue,
	'b=s@{1,}' => \@blast,
	'top=i' => \$top
);

## Parsing TCDB fasta
open TCDB, "<$tcdb";
my %tcdb;
while  (my $line = <TCDB>){
	chomp $line;
	if ($line =~ /^>(\S+)\s+(.*)$/){
		my $tag = $1; my $desc = $2;
		my $def = $tag; $def =~ s/^gnl\|TC-DB\|\S+\|//;
		$tcdb{$tag}[0] = $desc;
		$tcdb{$tag}[1] = $def;
	}
}

## Parsing blastp file
while (my $file = shift@blast){
	open IN, "<$file";
	open OUT, ">$file.parsed.top$top";
	print OUT "#locus\tE-value\tTCID family\tTCDB header\tTCDB description\n";
	my %seen = ();
	while (my $line = <IN>){
		chomp $line;
		my @cols = split("\t", $line);
		if ($cols[10] <= $evalue){
			if (exists $seen{$cols[0]}){$seen{$cols[0]} += 1;}
			else {$seen{$cols[0]} = 1;}
			if($seen{$cols[0]} <= $top){
				print OUT "$cols[0]\t$cols[10]\t$tcdb{$cols[1]}[1]\t$cols[1]\t$tcdb{$cols[1]}[0]\n";
			}
		}
	}
}