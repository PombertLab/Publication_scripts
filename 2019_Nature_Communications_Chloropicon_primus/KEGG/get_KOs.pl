#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = 'USAGE = get_KOs.pl -ko *.ko -o listname';
die "\n$usage\n\n" unless @ARGV;

my @ko;
my $output;

GetOptions(
	'ko=s@{1,}' => \@ko,
	'o=s' => \$output
);

open OUT, ">$output.ko";
while (my $file = shift@ko){
	open IN, "<$file";
	my %kos;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^\S+\tko:(K\d+)/){
			my $ko = $1;
			if (exists $kos{$ko}){next;}
			else{print OUT "$ko\n"; $kos{$ko} = $ko;}
		}
	}
}
