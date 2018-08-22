#!/usr/bin/perl

use strict;
use warnings;

my $usage = "USAGE = subsets.pl NUMBER_OF_SEQUENCES_PER_SUBSET FASTA_FILE\nEXAMPLE = subsets.pl 1000 CCMP1205.fasta\n";
die "\n$usage\n\n" unless @ARGV;

my $num = $ARGV[0];
my $file = $ARGV[1]; open FASTA, "<$file";

## Filling hash of sequences
my %sequences = (); my @seqs = ();
my $key = undef;
while (my $line = <FASTA>){
	chomp $line;
	if ($line =~ /^>(.*)$/){$key = $1; push(@seqs, $1);}
	else{$sequences{$key} .= $line;}
}
my $size = scalar(@seqs);
## Creating subsets
my $suffix = 0;
my $end = $size - $num;
for (my $i = 1; $i <= $end; $i += $num){
	$suffix++;
	open OUT, ">${file}_$suffix";
	for (1..$num){
		my $tmp = shift@seqs;
		print OUT ">$tmp\n";
		my @SEQ= unpack ("(A60)*", $sequences{$tmp});
		while (my $block = shift@SEQ){print OUT "$block\n";}
	}
}
$suffix++;
open OUT, ">${file}_$suffix";
while (my $tmp = shift@seqs){
	print OUT ">$tmp\n";
	my @SEQ= unpack ("(A60)*", $sequences{$tmp});
	while (my $block = shift@SEQ){print OUT "$block\n";}
}