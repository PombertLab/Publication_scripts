#!/usr/bin/perl

use strict;
use warnings;

my $usage = "USAGE = double_ko.pl ko00000.keg Verified_products_ALL.list CCMP_user_ko_definition.txt";
die "\n$usage\n\n" unless @ARGV;

open KEG, "<$ARGV[0]";
open CCMP, "<$ARGV[1]";
open KO, "<$ARGV[2]";

## Creating KEGG hash
my %kegg;
while (my $line = <KEG>){
	chomp $line;
	if ($line =~ /^D\s+(K\d+)\s+(.*)$/){
		my $ko = $1;
		my $def = $2;
		if (exists $kegg{$ko}){next;}
		else {$kegg{$ko} = $def;}
	}
}

## Creating CCMP prodcut list
my %prod;
while (my $line = <CCMP>){
	chomp $line;
	if ($line =~ /^#/){next;}
	elsif ($line =~ /^(\w+)\s+(.*)$/){$prod{$1} = $2;}
}

## Working on GHOSTKOALA KO list
open ALL, ">all.txt";
while (my $line = <KO>){
	chomp $line;
	my @tmp = split("\t", $line);
	if ($tmp[1] ne '') {
		if ($tmp[4] =~ /^(K\d+)/){
			print ALL "$tmp[0]\t$prod{$tmp[0]}\tPrimary and secondary KOs found\t$tmp[1]\t$kegg{$tmp[1]}\tKO #2\t$tmp[4]\t$kegg{$1}\n";
		}
		else{
			print ALL "$tmp[0]\t$prod{$tmp[0]}\tPrimary KO found\t$tmp[1]\t$kegg{$tmp[1]}\n";
		}
	}
	elsif ($tmp[4] =~ /^(K\d+)/){
		print ALL "$tmp[0]\t$prod{$tmp[0]}\tSecondary KO found\t$tmp[4]\t$kegg{$1}\n";
	}
	else{
		print ALL "$tmp[0]\t$prod{$tmp[0]}\tNo KO found\n"
	}
}