#!/usr/bin/perl
## Pombert Lab, 2022
my $name = 'methyldib.pl';
my $version = '0.1a';
my $updated = '2023-04-06';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
UPDATED		$updated
SYNOPSIS	Generates methyl distributions from Megalodon BED files for
        easy plotting with MS Exel and Circos.

COMMAND		$name \\
		  -f ref.fasta \\
		  -b *.bed \\
		  -o output_directory \\
		  -w 1000 \\
		  -s 500 \\
		  -c

-f (--fasta)	Reference fasta file
-b (--bed)      BED file to process
-o (--outdir)	Output directory [Default: ntBiases]
-w (--winsize)	Sliding window size [Default: 1000]
-s (--step)	Sliding window step [Default: 500]
-c (--circos)	Output files for Circos plotting
OPTIONS
die "\n$usage\n" unless @ARGV;

my $fasta;
my $bed;
my $outdir = 'ntBiases';
my $winsize = 1000;
my $step = 500;
my $circos;
GetOptions(
	'f|fasta=s' => \$fasta,
    'b|bed=s' => \$bed,
	'o|outdir=s' => \$outdir,
	'w|winsize=i' => \$winsize,
	's|step=i' => \$step,
	'c|circos' => \$circos
);

### Check if output directory / subdirs can be created
$outdir =~ s/\/$//;
unless (-d $outdir) {
	make_path($outdir,{mode => 0755})  or die "Can't create $outdir: $!\n";
}

### Working on reference FASTA file
my %sequences;
my $seq;

open FASTA, "<", $fasta or die "Can't open $fasta: $!\n";
while (my $line = <FASTA>){
    chomp $line;
    if ($line =~ /^>(\S+)/){
        $seq = $1;
    }
    else {
        $sequences{$seq} .= $line;
    }
}
close FASTA;

### Creating database of loci
my %loci_db;
foreach my $key (sort (keys %sequences)){
    my $sequence = $sequences{$key};
    my $len = length $sequence;
    for my $base (1..$len){
        $loci_db{$key}{$base} = 0;
    }
}

### Working on BED file
open BED, "<", $bed or die "Can't open $bed: $!\n";
while (my $line = <BED>){
    chomp $line;
    my @data = split("\t", $line);
    my $contig = $data[0];
    my $locus = $data[1];
    my $base_mod_per = $data[-1];
    $loci_db{$contig}{$locus} = $base_mod_per;
}

### Creating distribution
my ($basename) = fileparse($fasta);
my ($fileprefix) = $basename =~ /^(\w+)/;
my $circosfile = $outdir.'/'.$fileprefix.'.circos';

if ($circos){
	open CIRCOS, ">", $circosfile or die "Can't create $circosfile: $!\n";
	print CIRCOS "#chr START END  % methyl\n";
}

foreach my $contig (sort (keys %loci_db)){

    my $outfile = $outdir.'/'.$fileprefix.'.'.$contig.'.tsv';
	open METH, ">", $outfile or die "Can't create $outfile: $!\n";

    ### Sliding windows
    my $len = scalar (keys %{$loci_db{$contig}});

    my $x;
    my $sum;

    for ($x = 1; $x <= ($len - $winsize); $x += $step){

        for my $num ($x..($winsize+$x)){
            $sum += $loci_db{$contig}{$num};
        }

        my $start = $x - 1;
        my $end = $x + $winsize - 1;
		my $average = ($sum/$winsize);
        $sum = 0;
        print METH "$start\t$average\n";
		
		if ($circos){
			print CIRCOS "$contig $start $end $average\n";
		}

    }

	## left_over
	my $count = 0;
	for my $num ($x..$len){
		$count++;
		$sum += $loci_db{$contig}{$num};
	}

	my $start = $x - 1;
	my $end = $x + $count - 2;
	my $average = ($sum/$count);
    $sum = 0;
    print METH "$x\t$average\n";

	if ($circos){
		print CIRCOS "$contig $start $end $average\n";
	}

    close METH;

}
