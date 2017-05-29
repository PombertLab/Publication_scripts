#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2017

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $options = <<'OPTIONS';

EXAMPLE = cat_phy.pl -otu 26 -phy *.phy -o supermatrix.phy

-otu	# Number of expected OTUs [NOTE: this script does not handle missing data]
-phy	# Phylip files to be concatenated
-o		# Supermatix output file name

OPTIONS
die $options unless @ARGV;

## Defining variables
my $OTU;
my @phy;
my $out;
GetOptions(
	'otu=i' => \$OTU,
	'phy=s@{1,}' => \@phy,
	'out=s' => \$out
);

## Creating hash of concatenated sequences
my %py= ();
my $key = undef;
my $seq = undef;
while (my $phy = shift@phy){
	open IN, "<$phy" or die "Can\'t open $phy\n";
	my $taxa = '0';
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /\s+(\d+)\s+\d+/){$taxa = $1;}
		if ($taxa == $OTU){
			if ($line =~ /(\S+?)\@\S+\s+(.+)/){
				$key = $1;
				$seq = $2;
				$py{$key}.=$seq;
			}
		}
	}
}
### End of hash
my $taxa = scalar(keys(%py));	### number of OTUs
my $len = length$py{$key};		### length of concatenated matrix

open PHY, ">$out" or die "Can\'t create file $phy; please check file permissions\n";
print PHY " $taxa $len\n";
for (keys %py){print PHY "$_  $py{$_}\n";}

open FASTA, ">$out" or die "Can\'t create file $out; please check file permissions\n";
for (keys %py){print FASTA ">$_\n$py{$_}\n";}

