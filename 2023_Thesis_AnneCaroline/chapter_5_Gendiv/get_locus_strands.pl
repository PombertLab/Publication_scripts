#!/usr/bin/env perl
## Pombert Lab, Illinois Tech, 2023

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";
$0 -gbk *.gbk -gff *.gff -o locus_tag_list.txt
USAGE
die "\n$usage\n" unless @ARGV;

my @GBK;
my @GFF;
my $outfile = 'locus_tag_list.txt';
GetOptions(
	'gbk=s@{1,}' => \@GBK,
	'gff=s@{1,}' => \@GFF,
	'o|out=s' => \$outfile
);

my %loci;

while (my $gbk = shift @GBK){

	open GBK, '<', $gbk or die $!;

	my $gene_strand;
	my $locus_tag;

	while (my $line = <GBK>){

		chomp $line;

		if ($line =~ /^     gene            complement/){
			$gene_strand = '-';
		}
		elsif ($line =~ /^     gene            /){
			$gene_strand = '+';
		}
		if ($line =~ /^\s+\/locus_tag=\"(\w+)/){
			$locus_tag = $1;
		}

		$loci{$locus_tag} = $gene_strand;

	}

}

while (my $gff = shift @GFF){

	open GFF, '<', $gff or die $!;

	my $gene_strand;
	my $locus_tag;

	while (my $line = <GFF>){

		chomp $line;

		if ($line =~ /gene\t\d+\t\d+\t\.\t(\+|\-)\t\.\tID\=(\w+)/){
			$gene_strand = $1;
			$locus_tag = $2;
			$loci{$locus_tag} = $gene_strand;
		}	

	}

}


open OUT, '>', $outfile or die $!;

foreach my $locus (sort (keys %loci) ){

	print OUT "$locus\t$loci{$locus}\n";

}