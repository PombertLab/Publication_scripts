#!usr/bin/perl

open(IN, $ARGV[0]);
open OUT, ">volvox.fasta";

while(my $line = <IN>){
	chomp $line;
	if ($line =~ /^>jgi\|Volca1\|(\d{1,})/){
		my $gene=$1;
		print OUT ">VOLCADRAFT_$gene\n";
	}else{
		print OUT "$line\n";
	}
}