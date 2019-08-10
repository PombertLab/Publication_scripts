#!/usr/bin/perl

## Pombert Lab, IIT 2014
## Reorders the locus_tags in proper order in the EMBL files.

use strict;
use warnings;

system "mkdir GOOD";

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.embl//;
	open OUT, ">$file.ordered";
	open SEQ, ">$file.seq";
	$file =~ s/contig-//;
	
	my %gene = ();
	my %CDS = ();
	my @todolist = ();
	
	sub numSort {if ($a < $b) { return -1; }elsif ($a == $b) { return 0;}elsif ($a > $b) { return 1; }}
	
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /FT\s+gene\s+(\d+)/){
			my $start = $1;
			push (@todolist, $start);
			$gene{$start}=$line;
		}
		elsif ($line =~ /FT\s+gene\s+complement\((\d+)/){
			my $start = $1;
			push (@todolist, $start);
			$gene{$start}=$line;
		}
		elsif ($line =~ /FT\s+CDS\s+(.*)/){
			my @array = split(',',$1);
			my @exons= ();
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $start = $1;
					my $stop = $2;
					push (@exons, $start);
					push (@exons, $stop);
				}
			}
			my @ok = sort numSort @exons;			
			my $start = $ok[0];
			$CDS{$start}=$line;
		}
		elsif ($line =~ /FT\s+CDS\s+join\((.*)/){
			my @array = split(',',$1);
			my @exons= ();
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $start = $1;
					my $stop = $2;
					push (@exons, $start);
					push (@exons, $stop);
				}
			}
			my @ok = sort numSort @exons;			
			my $start = $ok[0];
			$CDS{$start}=$line;
		}
		elsif ($line =~ /FT\s+CDS\s+complement\((.*)/){
			my @array = split(',',$1);
			my @exons= ();
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $start = $1;
					my $stop = $2;
					push (@exons, $start);
					push (@exons, $stop);
				}
			}
			my @ok = sort numSort @exons;			
			my $start = $ok[0];
			$CDS{$start}=$line;
		}
		elsif ($line =~ /FT\s+CDS\s+complement\(join\((.*)/){
			my @array = split(',',$1);
			my @exons= ();
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $start = $1;
					my $stop = $2;
					push (@exons, $start);
					push (@exons, $stop);
				}
			}
			my @ok = sort numSort @exons;			
			my $start = $ok[0];
			$CDS{$start}=$line;
		}
		else{print SEQ "$line\n";}
	}
	
	my @todo = sort numSort @todolist;
	my $count = 10;
	
	while (my $feat = shift@todo){
		chomp $feat;
		print OUT "$gene{$feat}\n";
		print OUT "FT                   /locus_tag=\"DI09_$file"."p$count\"\n";
		print OUT "$CDS{$feat}\n";
		$count+=10;
	}
	system "cat contig-$file.ordered contig-$file.seq > GOOD/contig-$file.embl";
	system "rm *.ordered *.seq";
}