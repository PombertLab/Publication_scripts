#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my @fasta;
my @loci;
my $outdir = './';
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'l|loci=s@{1,}' => \@loci,
	'o|out=s' => \$outdir
);

unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
}

########################################

my %sequences;

while (my $fasta = shift @fasta){

	open FASTA, "<", $fasta or die $!;

	my $locus;

	while (my $line = <FASTA>){

		chomp $line;

		if ($line =~ /^>(\S+)/){
			$locus = $1;
			$sequences{$locus}{'header'} = $line;
		}
		else {
			$sequences{$locus}{'sequence'} .= $line;
		}

	}

}

########################################

while (my $coordinates = shift @loci){

	## struct = CP119067.1:51565-52986
	my ($locus,$coord) = split (":", $coordinates);
	my ($start,$end) = split ("-", $coord);
	my $len = $end - $start + 1; 

	my $subset = substr($sequences{$locus}{'sequence'}, $start - 1, $len);

	my $outfile = $outdir.'/'.$coordinates.'.fasta';
	open my $fh, '>', $outfile or die;

	print $fh $sequences{$locus}{'header'}."\n";

	my @seq = unpack ("(A60)*", $subset);
		while (my $seq = shift@seq){
			print $fh $seq."\n";
	}

}