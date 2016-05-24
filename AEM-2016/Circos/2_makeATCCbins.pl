#!/usr/bin/perl

use strict;
use warnings;

open IN1, "<ATCC.vcf.contig1";
open IN2, "<ATCC.vcf.contig2";
open IN3, "<ATCC.vcf.contig3";
open IN4, "<ATCC.vcf.contig4";

open OUT1, ">ATCC.contig1.bin";
open OUT2, ">ATCC.contig2.bin";
open OUT3, ">ATCC.contig3.bin";
open OUT4, ">ATCC.contig4.bin";

my %hash1=();
my %hash2=();
my %hash3=();
my %hash4=();

while (my $h1 = <IN1>){
	chomp $h1;
	if ($h1 =~ /^gi\|\d+\|\w+\|(\w+)\S+\|\s(\d+)/){$hash1{$2}=$2;}
}
for my $pos1 (1..4003){
	if (exists $hash1{$pos1}){print OUT1 "1";}
	else {print OUT1 "0";}
}

while (my $h2 = <IN2>){
	chomp $h2;
	if ($h2 =~ /^gi\|\d+\|\w+\|(\w+)\S+\|\s(\d+)/){$hash2{$2}=$2;}
}
for my $pos2 (1..527678){
	if (exists $hash2{$pos2}){print OUT2 "1";}
	else {print OUT2 "0";}
}

while (my $h3 = <IN3>){
	chomp $h3;
	if ($h3 =~ /^gi\|\d+\|\w+\|(\w+)\S+\|\s(\d+)/){$hash3{$2}=$2;}
}
for my $pos3 (1..1391602){
	if (exists $hash3{$pos3}){print OUT3 "1";}
	else {print OUT3 "0";}
}

while (my $h4 = <IN4>){
	chomp $h4;
	if ($h4 =~ /^gi\|\d+\|\w+\|(\w+)\S+\|\s(\d+)/){$hash4{$2}=$2;}
}
for my $pos4(1..2178842){
	if (exists $hash4{$pos4}){print OUT4 "1";}
	else {print OUT4 "0";}
}