#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = "USAGE = ./clusters.pl -t TAXA -f FASTA -g gap_size\nEXAMPLE: ./clusters -t CCMP -f ../BLAST/FASTA/CCMP1205_proteins.fasta -g 50";
die "\n$usage\n\n" unless @ARGV;

my @list = ('AUXENO','BATHY','CCMP','CHLAMY','CHLORELLA','COCCO','GONIUM','HELICO','KLEBS','MCOMMODA','MPUSI','OLUCI','OTAURI','VOLVOX');
my $taxa;
my $fasta;
my $gap = 50;

GetOptions(
	't|taxa=s' => \$taxa,
	'f|fasta=s' => \$fasta,
	'g|gap=i' => \$gap
);
while (my $list = shift@list){
	system "../get_synteny.pl -q ../LISTS/$taxa.list -b ../BLAST/${taxa}vs/${taxa}vs$list.blastp.6 -s ../LISTS/$list.list -o ${taxa}.vs.$list -g $gap";
}
my $prot_num = `grep -c '>' $fasta`;
chomp $prot_num;
open OUT, ">$taxa.gene_pairs.gap$gap.tsv";
print OUT "#Species_compared\tGene_pairs\t$taxa protein total = $prot_num\n";
system "grep -c -v '#' $taxa.*.gap$gap.clusters | sed 's/:/\t/g'";
system "grep -c -v '#' $taxa.*.gap$gap.clusters | sed 's/:/\t/g' >> $taxa.gene_pairs.gap$gap.tsv";
