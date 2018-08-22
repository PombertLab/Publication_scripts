#!/usr/bin/perl
## Pombert Lab, 2017

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = "USAGE = GC_content_to_Circos.pl -f *.fasta";
die "\n$usage\t# Type -h for help\n\n" unless @ARGV;

my $options = <<'OPTIONS';

EXAMPLE: GC_content_to_Circos.pl -f *.fasta -s 500 -w 1000

-h (-help)	Displays list of options	
-f (--fasta)	Input files in fasta format
-o (--ouput)	Output file names prefix [Default: genome]
-c (--color)	Color for genotype [Default: black]
-s (--step)	Size of the steps between windows [Default: 500]
-w (--window)	Width of the sliding windows [Default: 1000]

OPTIONS

my $help;
my @fasta;
my $output = 'genome';
my $color = 'black';
my $slide = 500; my $window = 1000;

GetOptions(
	'h|help' => \$help,
	'f|fasta=s@{1,}' => \@fasta,
	'o|ouput=s' => \$output,
	'c|color=s' => \$color,
	's|slide=i' => \$slide,
	'w|window=i' => \$window
);

die $options if $help;

my %sequences;
my @seq_list;
my $seq;

## Populating list of sequences
while (my $fasta = shift @fasta) {
	open IN, "<$fasta" or die "Can\'t open $fasta";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\w+).*$/){$seq = $1; push (@seq_list, $seq);}
		else {$sequences{$seq} .= $line;}
	}
	close IN;
}

## Creating a "karyotype" file for Circos
open GC, ">$output.gc"; print GC '#chr START END GC_PERCENTAGE'."\n";
open KAR, ">$output.genotype"; print KAR '#chr - ID LABEL START END COLOR'."\n";
my $id = 0; my $mmax; my $mmin;

## Working on sequences
while (my $sequence = shift @seq_list){
	my $len = length($sequences{$sequence});
	my $terminus = $len - 1;
	my $max; my $min;
	$id ++;
	print KAR "chr - $sequence $id 0 $terminus $color\n";
	print "Working on $sequence ... ";
	for(my $i = 0; $i <= $len - ($window-1); $i += $slide) {
		my $sub = substr($sequences{$sequence}, $i, $window);
		#print "$sub\n"; ## Debugging
		my $GC = 0;
		my @nt = unpack("(A1)*", $sub);
		my $size = scalar(@nt);
		while (my $base = shift@nt){if ($base =~ /[GCgc]/){$GC++;}}
		my $percent = sprintf("%.1f", ($GC/$size)*100);
		my $end = $i + $size - 1;
		print GC "$sequence $i $end $percent\n";
		if (!defined $max){$max = $percent;}
		if (!defined $mmax){$mmax = $percent;}
		if (!defined $min){$min = $percent;}
		if (!defined $mmin){$mmin = $percent;}
		if ($percent > $max){$max = $percent;}
		if ($percent > $mmax){$mmax = $percent;}
		if ($percent < $min){$min = $percent;}
		if ($percent < $mmin){$mmin = $percent;}
	}
	print "Min/Max GC = $min / $max \%\n";
}
print "Min GC percent = $mmin \%\n";
print "Max GC percent = $mmax \%\n";
close GC;
close KAR;
exit;
