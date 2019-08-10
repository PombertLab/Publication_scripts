#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl partial_mRNAs_tags_to_tbl *.tbl';

die $usage unless @ARGV;



while (my $tbl= shift @ARGV) {

	open IN, "<$tbl" or die "cannot open $tbl";
	$tbl =~ s/\.tbl$//;
	open OUT, ">$tbl.tbl2";
	
	local $/;	## slurp the content of the full file into the scalar ##
	my $complete_file = <IN>; ## slurp the content of the full file into the scalar ##

	$complete_file =~ s/\t(\d+)(\n\t\t\ttranscript_id)/\t>$1$2/mg;
	$complete_file =~s/\t(\d+)(\tmRNA\n\t\t\ttranscript_id)/\t>$1$2/mg;
	
	print OUT "$complete_file";
	
	close IN;
	
}
