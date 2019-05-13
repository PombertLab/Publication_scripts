#!usr/bin/perl

open(IN, $ARGV[0]);
open OUT, ">chlamy.fasta";

while(my $line = <IN>){
	chomp $line;
	if ($line =~ /^>jgi\|Chlre3\|(\d{1,})/){
		my $gene=$1;
		print OUT ">CHLREDRAFT_$gene\n";
	}else{
		print OUT "$line\n";
	}
}