#!/usr/bin/perl
## Pombert Lab, 2018
## Quick script to generate tables from the SNAP summary files

use strict;
use warnings;

die "USAGE = cat_summaries.pl products.txt summary.*\n" unless @ARGV;

## Products database
my $products = shift@ARGV;
open PD, "<$products";
my %prod;
while (my $line = <PD>){
	chomp $line;
	if ($line =~ /^(\S+)\t(.*)$/){$prod{$1} = $2;}
}

my %species = (
	'AEWD_' => 'Encephalitozoon_cuniculi_EC1', 'AEWQ_' => 'Encephalitozoon_cuniculi_EC2', 'AEWR_' => 'Encephalitozoon_cuniculi_EC3',
	'ECU01_' => 'Encephalitozoon_cuniculi_GB-M1',	'ECU02_' => 'Encephalitozoon_cuniculi_GB-M1', 'ECU03_' => 'Encephalitozoon_cuniculi_GB-M1', 'ECU04_' => 'Encephalitozoon_cuniculi_GB-M1', 'ECU05_' => 'Encephalitozoon_cuniculi_GB-M1',
	'ECU06_' => 'Encephalitozoon_cuniculi_GB-M1', 'ECU07_' => 'Encephalitozoon_cuniculi_GB-M1', 'ECU08_' => 'Encephalitozoon_cuniculi_GB-M1', 'ECU09_' => 'Encephalitozoon_cuniculi_GB-M1', 'ECU10_' => 'Encephalitozoon_cuniculi_GB-M1',
	'ECU11_' => 'Encephalitozoon_cuniculi_GB-M1',
	'EHEL_' => 'Encephalitozoon_hellem_ATCC_50504', 'KMI_' => 'Encephalitozoon_hellem_Swiss', 'Eint_' => 'Encephalitozoon_intestinalis_ATCC_50506', 'EROM_' => 'Encephalitozoon_romaleae_SJ-2008',
	'CWI42_' => 'Ordospora_colligata_FI-SK-17-1', 'CWI41_' => 'Ordospora_colligata_GB-EP-1', 'CWI40_' => 'Ordospora_colligata_NO-V-7', 'M896_' => 'Ordospora_colligata__OC4',
	'CWI36_' => 'Hamiltosporidium_magnivora_BE-OM-2', 'CWI37_' => 'Hamiltosporidium_tvaerminnensis_FI-OER-3-3',  'CWI38_' => 'Hamiltosporidium_tvaerminnensis_IL-G-3', 'CWI39_' => 'Hamiltosporidium_magnivora_IL-BN-2',
	'Napis_' => 'Nosema_apis_BRL01', 'NBO_' => 'Nosema_bombycis_CQ1', 'NCER_' => 'Nosema_ceranae_BRL01'	
);

## Working on summaries
open OUT, ">summaries_table.tsv";
print OUT "Orthogroup\tSpecies_1\tseq1\tProduct(seq1)\tSpecies_2\tseq2\tProduct(seq2)\tds\tdn\tds\/dn\tdn\/ds\n";
while (my $summary = shift@ARGV){
	open SUM, "<$summary";
	my $og = $summary; $og =~ s/^summary.//; $og =~ s/.NT$//; $og =~ s/.AA$//;
	while (my $line = <SUM>){
		chomp $line;
		if ($line =~ /^Compare/){next;}
		elsif ($line =~ /^Averages/){next;}
		else{
			my @col = split(/\s+/, $line);
			# $col[2] = sequence 1
			# $col[3] = sequence 2
			# $col[10] = ds; $col[11] = dn; $col[12] = ds/dn;
			my $sp1; my $sp2;
			if ($col[2] =~ /^(\S+_)/){$sp1 = $1;}
			if ($col[3] =~ /^(\S+_)/){$sp2 = $1;}
			my $dsdn;
			if (($col[10] eq 'NA') || ($col[11] eq 'NA')){$dsdn = 'NA';}
			elsif (($col[10] <= 0) || ($col[11] <= 0)){$dsdn = 'NA';}
			elsif (($col[10] > 0) && ($col[11] > 0)){$dsdn = ($col[11]/$col[10]); $dsdn = sprintf("%.4f", $dsdn);}
			print OUT "$og\t$species{$sp1}\t$col[2]\t$prod{$col[2]}\t$species{$sp2}\t$col[3]\t$prod{$col[3]}\t$col[10]\t$col[11]\t$col[12]\t$dsdn\n";
		}
	}
}