#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl script KEGG_file Table_KO_list';

open IN1, "<$ARGV[0]";
open IN2, "<$ARGV[1]";
open OUT, ">file.KOs";
	
my %orthology = ();
	
while (my $line = <IN1>){
	chomp $line;
	if ($line =~ /^ota:(Ot\d{2}g\d{5}).*(K\d{5})/){
		my $protein = $1;
		my $ko = $2;
		push (@{$orthology{$ko}}, $protein);
	}
}
while (my $kos = <IN2>){
	chomp $kos;
	if ($kos =~ /^\[KO:(K\d{5})\]/){
		my $target = $1;
			if (scalar(@{$orthology{$target}}) == 1){
				print OUT "$target\t"."$orthology{$target}[0]\n";
			}
			elsif (scalar(@{$orthology{$target}}) >= 2){
				my $nums = scalar(@{$orthology{$target}});
				my $end = $nums - 1;
				print OUT "$target\t";
				foreach my $count (0..$end-1){
					print OUT "$orthology{$target}[$count], ";
				}
				print OUT "$orthology{$target}[$end]\n";
			}
	}
}