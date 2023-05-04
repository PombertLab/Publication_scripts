#!/usr/bin/perl
## Pombert Lab, Illinois Tech, ACMdS, 2022

use strict;
use warnings;
use Chart::Gnuplot;
use Getopt::Long qw(GetOptions);
use File::Path qw(mkpath);

my $name = 'gc_plot.pl';
my $version = 0.3;
my $update = "2022-01-18";

my $options= <<"OPTIONS";
NAME    $name
VERSION $version
UPDATE  $update

OPTIONS:
-fa (--fasta)   (Multi)fasta file
-w (--winsize)  Size of window [Default:1000]
-s (--species)  Species name to put in plot title: e.g. 'Encephalitozoon hellem'
-f (--format)   Plot format [Default: eps]
-o (--outdir)   Output directory name [Default: fig]

OPTIONS
die "\n$options\n" unless @ARGV;

my $fasta;
my $outdir = 'fig';
my $winsize = 1000;
my $slide = int($winsize/2);
my $species = '';
my $format = 'eps';

GetOptions(

    'fa|fasta=s' => \$fasta,
    'w|winsize=i' => \$winsize,
    'o|outdir=s' => \$outdir,
    's|species=s' => \$species,
    'f|format=s' => \$format
);

## populate a hash of sequences
open FA, "<", "$fasta";

my %sequences;
my $seq_id;

while (my $line = <FA>){

    chomp $line;
    if ($line =~ /^>(\S+)/){
        $seq_id = $1;
    }
    else{
        $sequences{$seq_id} .= $line;
    }
}


## Check output directory
unless (-d $outdir){
    mkpath($outdir);
}

## Iterate through all sequences
foreach my $key (%sequences){

    ## GC content

    my $seq = $sequences{$key};
    my $length = length($seq);
    my $position = 0;
    my @gc;
    my @pos;

    for (my $i = $position; $i <= ($length - $winsize);  $i += $slide){
        my $DNA = substr($seq, $i, $winsize);
        my $count_c = $DNA =~ tr/Cc/Cc/;
        my $count_g = $DNA =~ tr/Gg/Gg/;
        my $start = $i + 1;
        my $end = $start + $winsize;
        my $gc_count = (($count_c + $count_g)/$winsize)*100;
        push(@pos, $start);
        push(@gc, $gc_count);

    }

    ## Data plot

    my $figname = "${key}_w${winsize}.$format";

    my $last = $length - $winsize;
    my $chart = Chart::Gnuplot->new(
        output => "${outdir}/${figname}",
        title  => "$species: $key",
        xlabel => "Position",
        ylabel => "GC%",
        xrange => [0, "$last"],
        yrange => [0, 100],
    );

    my $dataSet = Chart::Gnuplot::DataSet->new(
        xdata => \@pos,
        ydata => \@gc,
        title => "GC content",
        style => "linespoints",
    );
    
    # Plot the data set on the chart
    $chart->plot2d($dataSet);
    
    open OUT, ">", "${outdir}/${key}_w${winsize}.tsv";
    while (my $pos = shift@pos){
        my $gc = shift@gc;
        print OUT "$pos\t$gc\n";
    }

}