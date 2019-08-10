#!usr/bin/perl

open(IN, $ARGV[0]);
open OUT, ">lucima.fasta";

while(my $line = <IN>){
	chomp $line;
	if ($line =~ /^>jgi\|Ost\d{1,}_\d{1,}\|(\d{1,})/){
		my $gene=$1;
		print OUT ">OSTLU_$gene\n";
	}else{
		print OUT "$line\n";
	}
}