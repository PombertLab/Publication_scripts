#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'summarize_og.pl';
my $version = '0.1';
my $updated = '2022-07-07';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";

COMMAND		${name} -t Orthogroups.tsv -o table.tsv

USAGE
die "\n$usage\n" unless @ARGV;

my $tsv;
my $output;
GetOptions(
	't|tsv=s' => \$tsv,
	'o|output=s' => \$output
);

open OG, "<", $tsv or die "Can't open $tsv: $!\n";
open OUT, ">", $output or die "Can't create $output: $!\n";

print OUT "# Orthogroup\tECUN (JA071)\tEHEL (GPU96)\tEINT (GPK93)\n";

while (my $line = <OG>){

	chomp $line;

	unless ($line =~ /^Orthogroup/){
		
		my @data = split ("\t", $line);
		my $og = $data[0];
		my $ecu = $data[1];
		my $hel = $data[2];
		my $int = $data[3];

		print OUT "$og";

		foreach my $otu ($ecu, $hel, $int){
			if ($otu !~ '\S'){
				print OUT "\t0";
			}
			else {
				my @tmp = split (',', $otu);
				my $size = scalar (@tmp);
				print OUT "\t$size";
			}
		}
		print OUT "\n";
	}
}