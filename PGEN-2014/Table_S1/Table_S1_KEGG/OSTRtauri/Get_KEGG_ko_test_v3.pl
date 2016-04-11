#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl script KEGG_file Table_KO_list';

open IN1, "<$ARGV[0]";
open IN2, "<$ARGV[1]";
open OUT, ">KEGG.hash";

my %orthology = ();
	
while (my $line = <IN1>){
	chomp $line;
	if ($line =~ /^ota:(Ot\d{2}g\d{5}).*(K\d{5})/){
		my $protein = $1;
		my $ko = $2;
		push (@{$orthology{$ko}}, $protein);
		print OUT "$ko\t$protein\n";
	}
	elsif ($line =~ /^ota:(Ostap\w{2}\d+).*(K\d{5})/){
		my $protein = $1;
		my $ko = $2;
		push (@{$orthology{$ko}}, $protein);
		print OUT "$ko\t$protein\n";
	}
}
open IN3, "<KEGG.hash";
my %hash = ();
while (my $stuff = <IN3>){
	chomp $stuff;
	if ($stuff =~/(K\d{5})\t(Ot\d{2}g\d{5})/){
		$hash{$1} = $2;
	}
	elsif ($stuff =~ /^ota:(Ostap\w{2}\d+).*(K\d{5})/){
			$hash{$1} = $2;
	}
}

open OUT2, ">file.KOs";

while (my $kos = <IN2>){
	chomp $kos;
	if ($kos =~ /^\[KO:(K\d{5})\]/){
		my $target = $1;
			if (exists $hash{$target}){
				if (scalar(@{$orthology{$target}}) == 1){
					print OUT2 "$target\t"."$orthology{$target}[0]\n";
				}
				elsif (scalar(@{$orthology{$target}}) >= 2){
					my $nums = scalar(@{$orthology{$target}});
					my $end = $nums - 1;
					print OUT2 "$target\t";
					foreach my $count (0..$end-1){
						print OUT2 "$orthology{$target}[$count], ";
					}
					print OUT2 "$orthology{$target}[$end]\n";
				}
			}
			else{
				print OUT2 "$target\tNOT_FOUND\n";
			}
	}
}