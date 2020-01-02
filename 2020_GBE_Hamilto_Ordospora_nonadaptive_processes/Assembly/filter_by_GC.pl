#!/usr/bin/perl
## Pombert Lab, 2017
my $name = 'filter_by_GC.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Parses contigs by their GC content
USAGE		filter_by_GC.pl -min 20 -max 60 -fa *.fasta -out GC_summary.txt

OPTIONS:
-min		Desired minimum GC percentage [default: 20]
-max		Desired maximum GC percentage [default: 60]
-fa		Fasta file(s) to parse
-out		GC content summary output [default: GC_summary.txt]
OPTIONS
die "$options\n" unless @ARGV;

my $min = 20;
my $max = 60;
my @fasta;
my $output;
GetOptions(
	'min=i' => \$min,
	'max=i' => \$max,
	'fa|fasta=s@{1,}' => \@fasta,
	'out=s' => \$output
);

open OUT, ">$output";
while (my $fasta = shift@fasta){
	open IN, "<$fasta";
	$fasta =~ s/\.w+//;
	open MIN, ">$fasta.${min}-.GC.fasta";
	open MED, ">$fasta.${min}_to_${max}.GC.fasta";
	open MAX, ">$fasta.${max}+.GC.fasta";
	my %sequences;
	my $contig;
	my @contigs;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(.*)$/){
			push (@contigs, $1);
			$contig = $1;
		}
		else{
			$sequences{$contig} .= $line;
		}
	}
	while (my $cg = shift@contigs){
		my $gc = 0;
		my $length = length($sequences{$cg});
		for (my $i = 0; $i <= $length - 1;  $i++){
			my $nt = substr($sequences{$cg}, $i, 1);
			if ($nt eq 'G' || $nt eq 'g' || $nt eq 'C' || $nt eq 'c'){$gc++;}
		}
		my $percent = ($gc/$length)*100;
		if ($percent < $min){
			print OUT "UNDER\t$cg\tlength = $length\tGC percent = $percent\n";
			print MIN ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MIN "$sq\n";}
		}
		elsif ($percent > $max){
			print OUT "OVER\t$cg\tlength = $length\tGC percent = $percent\n";
			print MAX ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MAX "$sq\n";}
		}
		else{
			print OUT "BETWEEN\t$cg\tlength = $length\tGC percent = $percent\n";
			print MED ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MED "$sq\n";}
		}
	}
}
close OUT;
open IN, "<$output";
my $under = 0; my $between = 0; my $over = 0;
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^UNDER.*length = (\d+)/){$under+=$1;}
	elsif ($line =~ /^BETWEEN.*length = (\d+)/){$between+=$1;}
	elsif ($line =~ /^OVER.*length = (\d+)/){$over+=$1;}
}

my $sum = $under + $between + $over;
print "\nUnder ${min}\% GC\t\t$under bp\n";
print "In-between ${min}-${max}\% GC\t$between bp\n";
print "Over ${max}\ %GC\t\t$over bp\n";
print "Total\t\t\t$sum bp\n\n";