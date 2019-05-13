#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl partial_CDS_tags_to_tbl *.ends';

die $usage unless @ARGV;

## The idea here is to replace in the TBL files the values that are partial by the same values preceded by the proper > or < sign
## i.e. >$1, <$1, >$2, <$2, etc. Note that the (\s) after the value in the s/// operator is a must otherwise it would match any number starting with the same value(s)

while (my $ends= shift @ARGV) {

	open IN, "<$ends" or die "cannot open $ends";
	$ends =~ s/\.ends$//;
	open TBL2, "<$ends.tbl2";
	open OUT, ">$ends.tbl3";

	my $complete_file = undef;
	
	{
	local $/;	## slurp the content of the full file into the scalar ##
	$complete_file = <TBL2>; ## slurp the content of the full file into the scalar ##
	}
	
	while (my $line = <IN>) {
		
		chomp $line;
		
		### skip complete lines ###
			
		if ($line =~ m/^([+-])\t(H632_c\d+\.\d+)\t(complete)\t\w+\t\w+\t(\d+)\t(\d+)/) {
			next;
		}
		
		### work on plus strand ###
		
		elsif ($line =~ m/^(\+)\t(H632_c\d+\.\d+)\t(5end_missing)\t\w+\t\w+\t(\d+)\t(\d+)/) {
			my $strand = $1;
			my $locus = $2;
			my $status = $3;
			my $start = $4;
			my $end = $5;
					
				$complete_file =~ s/^($start)(\s)/<$1$2/mg;
		}
		elsif ($line =~ m/^(\+)\t(H632_c\d+\.\d+)\t(3end_missing)\t\w+\t\w+\t(\d+)\t(\d+)/) {
			my $strand = $1;
			my $locus = $2;
			my $status = $3;
			my $start = $4;
			my $end = $5;
		
				$complete_file =~ s/\t($end)(\s)/\t>$1$2/mg;
		}
		elsif ($line =~ m/^(\+)\t(H632_c\d+\.\d+)\t(both_ends_missing)\t\w+\t\w+\t(\d+)\t(\d+)/) {
			my $strand = $1;
			my $locus = $2;
			my $status = $3;
			my $start = $4;
			my $end = $5;
		
				$complete_file =~ s/^($start)(\s)/<$1$2/mg;
				$complete_file =~ s/\t($end)(\s)/\t>$1$2/mg;
		}
		
		### work on minus strand ###
		
		elsif ($line =~ m/^(\-)\t(H632_c\d+\.\d+)\t(5end_missing)\t\w+\t\w+\t(\d+)\t(\d+)/) {
			my $strand = $1;
			my $locus = $2;
			my $status = $3;
			my $start = $4;
			my $end = $5;
				
				$complete_file =~ s/^($start)(\s)/<$1$2/mg;
		}
		elsif ($line =~ m/^(\-)\t(H632_c\d+\.\d+)\t(3end_missing)\t\w+\t\w+\t(\d+)\t(\d+)/) {
			my $strand = $1;
			my $locus = $2;
			my $status = $3;
			my $start = $4;
			my $end = $5;
				
				$complete_file =~ s/\t($end)(\s)/\t>$1$2/mg;
		}
		elsif ($line =~ m/^(\-)\t(H632_c\d+\.\d+)\t(both_ends_missing)\t\w+\t\w+\t(\d+)\t(\d+)/) {
			my $strand = $1;
			my $locus = $2;
			my $status = $3;
			my $start = $4;
			my $end = $5;
				
				$complete_file =~ s/^($start)(\s)/<$1$2/mg;
				$complete_file =~ s/\t($end)(\s)/\t>$1$2/mg;
		}
		
	}
	
	print OUT "$complete_file";
	
	close IN;
	close TBL2;
	
}
