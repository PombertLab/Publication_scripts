#!/usr/bin/perl
## Pombert Lab, 2019
my $version = '0.1';
my $name = 'Ns_to_Excel.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generates a distribution of sites masked with RepeatMasker (N) using sliding windows
EXAMPLE		Ns_to_Excel.pl -f *.fasta -s 2500 -w 5000

-f (--fasta)	Input files in fasta format
-s (--step)	Size of the steps between windows [Default: 500]
-w (--window)	Width of the sliding windows [Default: 1000]
OPTIONS
die "$options\n" unless @ARGV;

my @fasta;
my $slide = 500;
my $window = 1000;

GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	's|slide=i' => \$slide,
	'w|window=i' => \$window
);

my %sequences;
my @seq_list;
my $seq;

## Populating list of sequences
while (my $fasta = shift @fasta) {
	open IN, "<$fasta" or die "Can\'t open $fasta";
	$fasta =~ s/\.\w+$//;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\w+).*$/){$seq = $1; push (@seq_list, $seq);}
		else {$sequences{$seq} .= $line;}
	}
	close IN;
}

## Working on sequences
my $id = 0; my $mmax; my $mmin;
while (my $sequence = shift @seq_list){
	my $len = length($sequences{$sequence});
	my $terminus = $len - 1;
	my $max; my $min;
	$id ++;
	print "Working on $sequence ... ";
	open NN, ">$sequence.s$slide.w$window.Ns.tsv";
	for(my $i = 0; $i <= $len - ($window-1); $i += $slide) {
		my $sub = substr($sequences{$sequence}, $i, $window);
		my $GC = 0;
		my @nt = unpack("(A1)*", $sub);
		my $size = scalar(@nt);
		while (my $base = shift@nt){if ($base =~ /[Nn]/){$GC++;}}
		my $percent = sprintf("%.1f", ($GC/$size)*100);
		my $end = $i + $size - 1;
		my $pos = $i + 1;
		print NN "$sequence\t$pos\t$percent\n";
		if (!defined $max){$max = $percent;}
		if (!defined $mmax){$mmax = $percent;}
		if (!defined $min){$min = $percent;}
		if (!defined $mmin){$mmin = $percent;}
		if ($percent > $max){$max = $percent;}
		if ($percent > $mmax){$mmax = $percent;}
		if ($percent < $min){$min = $percent;}
		if ($percent < $mmin){$mmin = $percent;}
	}
	print "Min/Max Ns = $min / $max \%\n";
}
print "Min Ns percent = $mmin \%\n";
print "Max Ns percent = $mmax \%\n";
close NN;
exit;
