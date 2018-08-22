#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

die "\nUSAGE = sort_SNPs.pl -min 35 -max 65 -vcf *.vcf\n\n" unless @ARGV;

my $min = 35;
my $max = 65;
my @VCF;
GetOptions(
	'min=i' => \$min,
	'max=i' => \$max,
	'v|vcf=s@{1,}' => \@VCF
);

while (my $file = shift@VCF){
	open IN, "<$file"; $file =~ s/.vcf$//;
	open OUT1, ">$file.sorted.$min.$max.vcf";
	open OUT2, ">$file.sorted.$min.$max.tsv";
	open OUT3, ">$file.unsorted.$min.$max.tsv"; 	## File to graph with excel or R
	open OUT4, ">$file.gaussian.$min.$max.tsv"; 	## Data for Gaussian plot 
	my @percents; my %percents;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else{
			my @columns = split("\t", $line);
			my @stats = split(":", $columns[9]);
			$stats[6] =~ s/%//; $stats[6] = sprintf("%.03d", $stats[6]);
			if ($columns[3] eq 'N'){next;} ## Skipping masked nucleotides
			elsif (($stats[6] >= $min) && ($stats[6] <= $max)){
				push (@percents, $stats[6]);
				$percents{$stats[6]} += 1;
				print OUT1 "$line\n";
			} 
		}
	}
	my @unsort = @percents;
	@percents = sort@percents;
	@percents = reverse@percents;
	while (my $pc = shift@percents){print OUT2 "$pc\n";}
	while (my $un = shift@unsort){print OUT3 "$un\n";}
	for (keys %percents){print OUT4 "$_\t$percents{$_}\n";}
}