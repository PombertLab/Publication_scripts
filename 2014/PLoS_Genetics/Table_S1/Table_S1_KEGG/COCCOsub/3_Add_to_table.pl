#!/usr/bin/perl

my $usage = 'perl script table_list hash.formatted';

open IN, "<$ARGV[0]";
open IN2, "<$ARGV[1]";

open OUT, ">add_to_list.txt";

## fill hash

%hash = ();

while (my $stuff = <IN2>){
	chomp $stuff;
	if ($stuff =~ /^(K\d{5})\t(.*)$/){
		$hash{$1} = $2;
	}
}

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^\d{5}_/){
		print OUT "$line\n";
	}
	elsif ($line =~ /^\[KO:(K\d{5})\]/){
		my $ko = $1;
		print OUT "$ko\t$hash{$ko}\n";
	}
}
		