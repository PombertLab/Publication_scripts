#!/usr/bin/perl

use strict;
use warnings;

open IN, "<Table_S2_v3_extended.1e-10.txt";
open KL, ">Shared_with_Klebsormidium.1e-10.tsv"; open KLU, ">Shared_uniquely_with_Klebsormidium.1e-10.tsv";
print KL "\t\t\t\t\tKlebsormidium flaccidum\t\t\t\n"; print KLU "\t\t\t\t\tKlebsormidium flaccidum\t\t\t\n";
open PR, ">Shared_with_Prasinophytes.1e-10.tsv"; open PRU, ">Shared_uniquely_with_Prasinophytes.1e-10.tsv";
open CO, ">Shared_with_CoreChlorophytes.1e-10.tsv"; open COU, ">Shared_uniquely_with_CoreChlorophytes.1e-10.tsv";
open SOLO, ">Found_only_in_Chloropicon.1e-10.tsv";

my %klebs;
my %prasinos;
my %cores;
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^A3770/){
		my @hits = split("\t", $line);
		if (($hits[5] ne '---') || ($hits[7] ne '---')){$klebs{$hits[0]} = $hits[1];}
		if (($hits[10] ne '---') || ($hits[12] ne '---') || ($hits[15] ne '---') || ($hits[17] ne '---') || ($hits[20] ne '---') || ($hits[22] ne '---') || ($hits[25] ne '---') || ($hits[27] ne '---') || ($hits[30] ne '---') || ($hits[32] ne '---')){$prasinos{$hits[0]} = $hits[1] ;}
		if (($hits[35] ne '---') || ($hits[37] ne '---') || ($hits[40] ne '---') || ($hits[42] ne '---') || ($hits[45] ne '---') || ($hits[47] ne '---') || ($hits[50] ne '---') || ($hits[53] ne '---') || ($hits[55] ne '---') || ($hits[57] ne '---') || ($hits[60] ne '---') || ($hits[62] ne '---') || ($hits[65] ne '---') || ($hits[67] ne '---')){$cores{$hits[0]} = $hits[1] ;}
		if (exists $klebs{$hits[0]}){
			print KL "$hits[0]\t$hits[1]\t";
			for (my $i = 5;  $i <= 8; $i++)  {print KL "$hits[$i]\t";}
			print KL "\n";
		}
		if (exists $prasinos{$hits[0]}){
			print PR "$hits[0]\t$hits[1]\t";
			for (my $i = 10;  $i <= 33; $i++)  {print PR "$hits[$i]\t";}
			print PR "\n";
		}
		if (exists $cores{$hits[0]}){
			print CO "$hits[0]\t$hits[1]\t";
			for (my $i = 35;  $i <= 68; $i++)  {print CO "$hits[$i]\t";}
			print CO "\n";
		}
		if ((exists $klebs{$hits[0]}) && (!exists $prasinos{$hits[0]}) && (!exists $cores{$hits[0]})){
			print KLU "$hits[0]\t$hits[1]\t";
			for (my $i = 5;  $i <= 8; $i++)  {print KLU "$hits[$i]\t";}
			print KLU "\n";
		}
		if ((exists $prasinos{$hits[0]}) && (!exists $klebs{$hits[0]}) && (!exists $cores{$hits[0]})){
			print PRU"$hits[0]\t$hits[1]\t";
			for (my $i = 10;  $i <= 33; $i++)  {print PRU "$hits[$i]\t";}
			print PRU "\n";
		}
		if ((exists $cores{$hits[0]}) && (!exists $klebs{$hits[0]}) && (!exists $prasinos{$hits[0]})){
			print COU "$hits[0]\t$hits[1]\t";
			for (my $i = 35;  $i <= 68; $i++)  {print COU "$hits[$i]\t";}
			print COU "\n";
		}
		if ((!exists $cores{$hits[0]}) && (!exists $klebs{$hits[0]}) && (!exists $prasinos{$hits[0]})){
			print SOLO "$hits[0]\t$hits[1]\n";
		}
	}
}