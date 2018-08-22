#!/usr/bin/perl
## Pombert Lab, 2017

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = "USAGE = TBL_to_list.pl -t *.tbl -ko *.ko.txt -o list";
die "\n$usage\t# Type -h for help\n\n" unless @ARGV;
my $options = <<'OPTIONS';

EXAMPLE: TBL_to_list.pl -t *.tbl -o list

-h (--help)	Displays list of options	
-t (--tbl)		Input files in TBL format
-ko		KEGG orthologs found with GHOSTKOALA
-o (--ouput)	Output file name prefix [Default: genome]

OPTIONS
my $help;
my @tbl;
my @ko;
my $output = 'genome';
GetOptions(
	'h|help' => \$help,
	't|tbl=s@{1,}' => \@tbl,
	'ko=s@{1,}' => \@ko,
	'o|ouput=s' => \$output
);
die $options if $help;

## Creating KOs database
my %KO;
while (my $file = shift@ko){
	open KO, "<$file";
	while (my $line = <KO>){
		chomp $line;
		if ($line =~ /^(\w+)\t(K\d+)/){$KO{$1} = $2;}
	}		
}

open OUT, ">$output.list";
while (my $file = shift@tbl){
	open IN, "<$file";
	$file =~ s/.tbl$//; $file =~ s/^.*\///;
	my $type;
	my $start;
	my $end;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^<?(\d+)\t>?(\d+)\t(\w+)/){
			$start = $1;
			$end = $2;
			$type = $3;
		}
		if (($line =~ /^\t\t\tlocus_tag\t(\w+)/) && ($type eq 'gene')){
			my $locus =$1;
			if ($start < $end) {
				if (exists $KO{$locus}){print OUT "$locus\t$file\t$start\t$end\tplus\t$KO{$locus}\n";}
				else {print OUT "$locus\t$file\t$start\t$end\tplus\tNA\n";}
			}
			elsif ($start > $end) {
				if (exists $KO{$locus}){print OUT "$locus\t$file\t$end\t$start\tminus\t$KO{$locus}\n";}
				else{print OUT "$locus\t$file\t$end\t$start\tminus\tNA\n";}
			}
		}
	}
}