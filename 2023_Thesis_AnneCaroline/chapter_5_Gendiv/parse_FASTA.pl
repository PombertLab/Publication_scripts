#!/usr/bin/env perl
## Pombert Lab, Illinois Tech, 2023

my $name = 'parse_FASTA.pl';
my $version = '0.1';
my $updated = '2023-05-20';

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Parses fasta file(s) by minimum size

COMMAND		${name} -m 1000 -f *.fasta -o outdir

OPTIONS:
-m (--min)		Minimum size [Default: 1000]
-f (--fasta)	Fasta file(s) to parse
-o (--out)		Output directory [Default: ./]
USAGE
die "\n$usage\n" unless @ARGV;

my @fasta;
my $minsize = 1000;
my $outdir = './';
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'm|min=i' => \$minsize,
	'o|out=s' => \$outdir	
);

unless (-d $outdir){
	mkdir($outdir, 0755) or die $!;
}

while (my $fasta = shift@fasta){

	open FASTA, '<', $fasta or die $!;
	my $basename = fileparse($fasta);
	$basename =~ s/\.\w+$//;
	my $outfile = $outdir.'/'.$basename.'.'.$minsize.'.fasta';
	open OUT, '>', $outfile or die $!;

	my %sequences;
	my $sequence;

	while (my $line = <FASTA>){
		
		chomp $line;
		if ($line =~ /^>/){
			$sequence = $line;
		}
		else {
			$sequences{$sequence} .= $line;
		}
	
	}

	## Sorting sequences longest -> shortest
	foreach my $seq (sort { length($sequences{$b}) <=> length($sequences{$a}) } keys %sequences){
		
		my $len = length $sequences{$seq};
		
		if ($len >= $minsize){

			print OUT $seq."\n";
			my @data = unpack("(A60)*", $sequences{$seq});
			while (my $tmp = shift@data){
				print OUT "$tmp\n";
			}

		}
	}

}