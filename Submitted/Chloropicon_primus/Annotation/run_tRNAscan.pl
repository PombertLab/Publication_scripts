#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	system "/opt/tRNAscan-SE/tRNAscan-SE $file > $file.tRNAs";
}
