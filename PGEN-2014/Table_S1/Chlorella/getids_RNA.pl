#!usr/bin/perl

use strict;
use warnings;

open PROT, "<vsChlorella.RNA.parsed.hash";
my %prot = ();

while (my $line = <PROT>) {
	chomp $line;
		if ($line =~ /^(CHLREDRAFT_\d+)\t(CBPS\S+)/) {
		my $key = $1;
		my $value = $2;
		$prot{"$key"} = "$value";
		}
}

while (my $file = shift @ARGV) {
	open IN, "<$file";
	open OUT, ">$file.ids";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(CHLREDRAFT_\d+)/){
			my $id = $1;
			print OUT "$prot{$id}\n";
		}else{
			print OUT "$line\n";
		}
	}
}