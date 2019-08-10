#!usr/bin/perl

open(IN, $ARGV[0]);
open OUT, ">otauri.fasta";

while(my $line = <IN>){
	chomp $line;
	if ($line =~ /gene=(Ot\d{2}g\d{5})/){
		my $gene=$1;
		print OUT ">$gene\n";
	}else{
		print OUT "$line\n";
	}
}