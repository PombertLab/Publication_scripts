#!/usr/bin/perl

## Pombert Lab, IIT, 2014
## Calculates the number of introns in protein-coding genes from EMBL files

use strict;
use warnings;

my $usage = 'USAGE = EMBLtoPROT *.embl';
die $usage unless @ARGV;

my $introns = 0;
my $prots = 0;

while (my $file = shift@ARGV){
	open IN, "<$file";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /FT   CDS             (\d+)..(\d+)/){$prots++;} ## Forward, single exon, skip
		elsif ($line =~ /FT   CDS             complement\((\d+)..(\d+)\)/){$prots++;} ## Reverse, single exon
		elsif ($line =~ /FT   CDS             join\((.*)\)/){ ## Forward, multiple exon
			my @array = split(',',$1);
			my @start = ();
			my @stop = ();
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $strt = $1;
					my $stp = $2;
					push (@start, $strt);
					push (@stop, $stp);
				}
			}
			my $asize = scalar(@start);
			my $int = ($asize-1);
			$introns += $int;
			$prots++;
		}
		elsif ($line =~ /FT   CDS             complement\(join\((.*)/){ ## Reverse, mutiple exon
			my @array = split(',',$1);
			my @start = ();
			my @stop = ();
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $strt = $1;
					my $stp = $2;
					unshift (@start, $strt);
					unshift (@stop, $stp);
				}
			}
			my $asize = scalar(@start);
			my $int = ($asize -1);
			$introns += $int;
			$prots++;
		}
	}
}
my $av = ($introns/$prots);
my $rd_av = sprintf("%.2f", $av);
print "The total protein count is: $prots\n";
print "The total intron count is: $introns\n";
print "The average intron per predicted protein is: $rd_av\n";