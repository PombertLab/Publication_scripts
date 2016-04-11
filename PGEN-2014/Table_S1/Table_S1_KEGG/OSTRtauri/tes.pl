#!/usr/bin/perl

open IN1, "<$ARGV[0]";
open OUT, ">test.out";

while (my $line = <IN1>){
	chomp $line;
	if ($line =~ /^ota:(Ot\d{2}g\d{5}).*(K\d{5})/){
		my $protein = $1;
		my $ko = $2;
		print OUT "$ko\t$protein\n";
	}
}