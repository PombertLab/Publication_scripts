#!/usr/bin/perl
## Pombert Lab 2018


use warnings;
use Getopt::Long qw(GetOptions);

my $usage= 'usage = convert6_matrix.pl -l .list -t .tsv

OPTIONS:
-g (--genes) list of gene
-a (--accession) list of accession
-t (--table) tsv table
';

die "\n$usage\n" unless @ARGV;

my $list; ##21 genes
my $table; ##input tsv table
GetOptions(
		'g|gene=s' => \$genelist,
		't|table=s' => \$table,
		'a|accession=s' => \$accs
);

##create list of genes
my @genes;
my %acc;

open IN, "<$genelist";
open OUT, ">matrix.txt";
while (my $line = <IN>){
	chomp $line;
	push @genes, "$line";
	print OUT ",$line";
}
print OUT "\n";
##create list of accessions
open ACC, "<$accs";
while (my $line = <ACC>){
	chomp $line;
	push @accession, "$line";
}

## create database of table
open TAB, "<$table";
while (my $line = <TAB>){
	chomp $line;
	if ($line =~ /^(\S+)\t(\S+)/){
		my $locus = $1;
		my $num = $2;
		$acc{$num}{$locus}+=1;
	}
}
##create matrix
while($li= shift@accession){
	print OUT "$li";
	foreach my $loop_variable (@genes){
		if (exists $acc{$li}{$loop_variable}){print OUT ",$acc{$li}{$loop_variable}";}
		else{print OUT ",0";}
	}
	print OUT "\n";
}


