#!/usr/bin/perl

## Sliding windows script for SNPs. Use binary 0 and 1 for the absence or presence of a SNPs at corresponding position.

use strict;
use warnings;

my $usage = 'perl SNPs_slide.pl *.bin';

die $usage unless @ARGV;

while (my $bin = shift @ARGV) {

	open IN, "<$bin" or die "cannot open $bin";
	$bin =~ s/\.bin$//;
	open OUT, ">$bin.windows";

	my $file = <IN>;
	chomp $file;
	my $len = length($file);
	
	## Sliding windows parameters
	
	my $winsize = 1000;	## replace with the desired width of the sliding windows.
	my $slide = 100;		## replace with the desired span of the slide, e.g. slide every nucleotide or every 10, 100 and so forth.
	
	## End of sliding windows parameters
	
	my $SNPs = 0;

	for(my $i = 0; $i <= $len-($winsize-1); $i += $slide) {
		
		my $window = substr($file, $i, $i+($winsize-1));
		my @bin = split('', $window);

			for (0..($winsize-2)) {
			
				if ($bin[$_] =~ /0/) {
					next;
				}
				elsif ($bin[$_] =~ /1/) {
					$SNPs++;
				}
			}
		
			for (($winsize-1)..($winsize-1)) {
			
				if ($bin[$_] =~ /1/) {
					$SNPs++;
					print OUT "$SNPs\n";
					$SNPs = 0;
				}
				else {
					print OUT "$SNPs\n";
					$SNPs = 0;
				}
			}
	}
}

close IN;
close OUT;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)