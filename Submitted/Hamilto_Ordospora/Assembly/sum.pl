#!/usr/bin/perl

use strict;
use warnings;

open IN, "<$ARGV[0]";
my $sum = 0;
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^UNDER.*lenght = (\d+)/){$sum+=$1;}
}
print "\nTOTAL = $sum bp\n\n";