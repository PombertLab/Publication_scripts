#!/usr/bin/perl
## Pombert lab, 2019
my $version = '0.1';
my $name = 'sort_methylation.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generates a distribution of methylated sites from the PacBio SMRT-LINK GFF methylation file using sliding windows
USAGE		sort_methylation.pl pacbio_methylation_file.gff

EXAMPLE		sort_methylation.pl -g pacbio_methylation_file.gff -f *.fasta -m 30 -o output

OPTIONS:
-h (-help)	Displays list of options	
-g (--gff)	GFF methylation file from SMRT-LINK
-f (--fasta)	Input files in fasta format
-m (--min)	Minimum identification quality value (Qv) [Default: 30]
-o (--output)	Output [Default: methylation]
-s (--step)	Size of the steps between windows [Default: 5000]
-w (--window)	Width of the sliding windows [Default: 10000]
OPTIONS
die "$options\n" unless @ARGV;

my $help;
my $gff;
my @fasta;
my $min = 30;
my $output = 'methylation';
my $slide = 5000; my $window = 10000;
GetOptions(
	'h|help' => \$help,
	'g|gff=s' => \$gff,
	'f|fasta=s@{1,}' => \@fasta,
	'm|min=i' => \$min,
	'o|ouput=s' => \$output,
	's|slide=i' => \$slide,
	'w|window=i' => \$window
);
die $options if $help;

## Populating list of methylated sites
open GFF, "<$gff" or die "Error. Can't open GFF methylation file $gff\n";
my %methyls; my $tot_methyls = 0;
while (my $line = <GFF>){
	chomp $line;
	if ($line =~ /^#/){next;}
	elsif ($line =~ /modified_base/){next;}
	else{
		my @columns = split("\t", $line);
		my $contig = $columns[0];
		my $type = $columns[2];
		my $position = $columns[3];
		my ($quality) = $columns[8] =~ /identificationQv=(\d+)/;
		if ($quality >= $min){$methyls{$contig}{$position} = $type; $tot_methyls++;}
	}
}
close GFF;
print "Found a total of $tot_methyls methylated sites at Qv >= $min\n";

## Populating list of sequences
my %sequences; my @seq_list; my $seq;
while (my $fasta = shift @fasta) {
	open IN, "<$fasta" or die "Can\'t open $fasta\n";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+).*$/){$seq = $1; push (@seq_list, $seq);}
		else {$sequences{$seq} .= $line;}
	}
	close IN;
}

## Creating sliding windows
open OUT, ">$output.tsv"; print OUT "#chr\tposition\tMethylation_count\n";
my $mmax; my $mmin;
while (my $sequence = shift @seq_list){
	my $len = length($sequences{$sequence});
	my $terminus = $len - 1;
	my $max; my $min;
	print "Working on $sequence ... ";
	for(my $i = 0; $i <= $len - ($window-1); $i += $slide) {
		my $met = 0;
		for ($i..($i+$window-1)){
			if (exists $methyls{$sequence}{$_}){$met++;}
		}
		my $end = $i + $window - 1;
		my $pos = $i + 1;
		print OUT "$sequence\t$pos\t$met\n";
		if (!defined $max){$max = $met;}
		if (!defined $mmax){$mmax = $met;}
		if (!defined $min){$min = $met;}
		if (!defined $mmin){$mmin = $met;}
		if ($met > $max){$max = $met;}
		if ($met > $mmax){$mmax = $met;}
		if ($met < $min){$min = $met;}
		if ($met < $mmin){$mmin = $met;}
	}
	print "Min/Max methylated sites per window = $min / $max\n";
}
print "Min methylated sites per window = $mmin \n";
print "Max methylated sites per window = $mmax \n";
print "Total number of methylated sites = $tot_methyls\n";
close OUT;
exit;

