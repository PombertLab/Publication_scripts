#!/usr/bin/perl

open IN, "<$ARGV[0]";
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^\S+\t(\S+)/){
		my $chromo = $1;
		open OUT, ">>$chromo.list";
		print OUT "$line\n";
		close OUT;
	}
}
