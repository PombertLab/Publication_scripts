#!/usr/bin/perl
## Pombert Lab, 2017

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## Regex to use for GFF parsing. Can also be entered from the command line.
my $regex = '^(\S+).*CDS\t(\d+)\t(\d+)\t\.\t([+-])\t[012.].*Genbank:(\w+\.\d+)';
# $1 => contig
# $2 => start , $3 => end
# $4 => strand, $5 => gene or protein

## Defining options
my $usage = "USAGE = gff_to_synteny.pl -i file.gff -o file.list -r 'regex'";
die "\n$usage\n\n" unless @ARGV;

my $options = << 'OPTIONS';

USAGE (simple) = gff_to_synteny.pl -i file.gff -o file.list
USAGE (advanced) = gff_to_synteny.pl -i file.gff -o file.list -r 'regex'

-h (--help)		Displays this help message
-i (--input)		Input file in GFF format
-o (--output)		Creates a list file for get_synteny.pl
-r (--regex)		Regular expression to parse GFF -		## GFF are not quite standard; the regex to capture the info has to be tweaked accordingly.
			Can also be set at the top of the script.	## The regex and its required paramaters to capture are located at the top of this script.

## Regex
$1 => contig
$2 => start
$3 => end
$4 => strand
$5 => gene or protein

OPTIONS

my $help;
my $input;
my $output;

GetOptions(
	'h|help' => \$help,
	'i|input=s' => \$input,
	'o|output=s' => \$output,
	'r|regex=s' => \$regex
);
die $options if $help;

## Parsing GFF
my %contigs;
my %genes;
open IN, "<$input" or die "Can't open input file : $input\n";
open OUT, ">$output" or die "Can't write to output file : $output\n";
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /$regex/){
		my $contig = $1;
		my $start = $2; my $end = $3;
		my $strand = $4; my $gene =$5;
		if (!exists $genes{$gene}){
			$genes{$gene} = $gene;
			$start = sprintf("%10d", $start); $end = sprintf("%10d", $end); ## Preventing number ordering SNAFU
			$contigs{$contig}{$start}[0] = $contig;
			$contigs{$contig}{$start}[1] = $start;
			$contigs{$contig}{$start}[2] = $end;
			$contigs{$contig}{$start}[3] = $strand;
			$contigs{$contig}{$start}[4] = $gene;
		}
	}
}
foreach my $cg (sort keys %contigs) { ## We must reorder genes in case that they were entered out of order in the GFF files
	my $position = 0;
	foreach my $pos (sort keys %{ $contigs{$cg} }) {
		$position ++;
		print OUT "$contigs{$cg}{$pos}[4]\t"; ## Printing gene
		print OUT "$contigs{$cg}{$pos}[0]\t"; ## Printing contig
		print OUT "$contigs{$cg}{$pos}[1]\t"; ## Printing $start
		print OUT "$contigs{$cg}{$pos}[2]\t"; ## Printing $end
		print OUT "$contigs{$cg}{$pos}[3]\t"; ## Printing $strand
		print OUT "$position\n"; ## Printing $strand
	}
}
