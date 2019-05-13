#!/usr/bin/perl

## Sliding windows script for G+C content.
## Uses a single DNA string files as input.

use strict;
use warnings;

my $usage = 'GC_slide.pl *.string';

die $usage unless @ARGV;

while (my $string = shift @ARGV) {

	open IN, "<$string" or die "cannot open $string";
	$string =~ s/\.string$//;
	open OUT, ">$string.gc";

	my $file = <IN>;
	chomp $file;
	my $len = length($file);
	
	## Sliding windows parameters
	
	my $winsize = 1000;	## replace with the desired width of the sliding windows.
	my $slide = 100;		## replace with the desired span of the slide, e.g. slide every nucleotide or every 10, 100 and so forth...
	
	## End of sliding windows parameters
	
	my $GC = 0;

	for(my $i = 0; $i <= $len-($winsize-1); $i += $slide) {
		
		my $window = substr($file, $i, $i+($winsize-1));
		my @nt = split('', $window);

			for (0..($winsize-2)) {
			
				if ($nt[$_] =~ /[A|T]/i) {
					next;
				}
				elsif ($nt[$_] =~ /[G|C]/i) {
					$GC++;
				}
			}
		
			for (($winsize-1)..($winsize-1)) {
			
				if ($nt[$_] =~ /[G|C]/i) {
					$GC++;
					my $GCpercent = ($GC/10);
					print OUT "$GCpercent\n";
					$GC = 0;
				}
				else {
					my $GCpercent = ($GC/10);
					print OUT "$GCpercent\n";
					$GC = 0;
				}
			}
	}
}

close IN;
close OUT;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)