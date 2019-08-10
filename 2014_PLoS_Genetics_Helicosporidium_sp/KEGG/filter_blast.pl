#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open CHLAMY, ">$file.chlamy";
	open COCCO, ">$file.cocco";
	open CHLORE, ">$file.chlore";
	open MICPU, ">$file.micpu";
	open OSTRE, ">$file.ostre";
	open PHYPA, ">$file.phypa";
	open VOLCA, ">$file.volca";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /Chlre4/){
			print CHLAMY "$line\n";
		}
		elsif ($line =~ /Coc_C169_1/){
			print COCCO "$line\n";
		}
		elsif ($line =~ /ChlNC64A_1/){
			print CHLORE "$line\n";
		}
		elsif ($line =~ /MicpuN3/){
			print MICPU "$line\n";
		}
		elsif ($line =~ /Ostta4/){
			print OSTRE "$line\n";
		}
		elsif ($line =~ /Phypa1_1/){
			print PHYPA "$line\n";
		}
		elsif ($line =~ /Volca1/){
			print VOLCA "$line\n";
		}
	}
}