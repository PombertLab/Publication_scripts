#!/usr/bin/env perl
## Pombert Lab, 2024

my $name = 'count_biases.pl';
my $version = '0.2';
my $updated = '2024-02-06';

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);

my $usage = <<"USAGE";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Calculates nucleotides biases in FASTA sequences

COMMAND     $0 \\
              -f *.fasta \\
              -o table.tsv \\
              -s summary.tsv

OPTIONS:
-f (--fasta)    Input FASTA file
-o (--out)      TAB-delimited output file (per sequence)
-s (--sum)      TAB-delimited output file (per file summary)
USAGE
die "\n$usage\n" unless @ARGV;

my @fasta;
my $outfile;
my $summary;
GetOptions(
    'f|fasta=s@{1,}' => \@fasta,
    'o|out=s' => \$outfile,
    's|sum=s' => \$summary
);

####### Outfile
open TSV, '>', $outfile or die "Can't create $outfile: $!\n";
open SUM, '>', $summary or die "Can't create $summary: $!\n";

print TSV 'File'."\t".'Chr'."\t".'%GC'."\t".'%AT'."\t";
print TSV '%GA'."\t".'%CT'."\t".'%GT'."\t".'%AC'."\t";
print TSV '%NN'."\t".'Total length'."\t".'Total NNs'."\n";

print SUM 'File'."\t".'%GC'."\t".'%AT'."\t";
print SUM '%GA'."\t".'%CT'."\t".'%GT'."\t".'%AC'."\t";
print SUM '%NN'."\t".'Total length'."\t".'Total NNs'."\n";

####### Parsing fasta files 
while (my $fasta = shift @fasta){

    open FASTA, '<', $fasta or die "Can't open $fasta: $!\n";
    my ($basename,$path) = fileparse($fasta);
    my %data;
    my $seqname;

    while (my $line = <FASTA>){
        chomp $line;
        if ($line =~ /^>(\S+)/){
            $seqname = $1;
        }
        else {
            $data{$seqname} .= $line;
        }
    }

    ## Metrics (total)
    my $tot_len;
    my $tot_gc;
    my $tot_at;
    my $tot_ga;
    my $tot_ct;
    my $tot_gt;
    my $tot_ac;
    my $tot_nn;

    foreach my $key (sort(keys(%data))){

        my $seq = $data{$key};
        my $len = length($seq);
        $tot_len += $len;

        ## Metrics (per contig/chromosome)
        my $gc = $seq =~ tr/GgCc//;
        my $at = $seq =~ tr/AaTt//;
        my $ga = $seq =~ tr/GgAa//;
        my $ct = $seq =~ tr/CcTt//;
        my $gt = $seq =~ tr/GgTt//;
        my $ac = $seq =~ tr/AaCc//;
        my $nn = $seq =~ tr/Nn//;

        $tot_gc += $gc;
        $tot_at += $at;
        $tot_ga += $ga;
        $tot_ct += $ct;
        $tot_gt += $gt;
        $tot_ac += $ac;
        $tot_nn += $nn;

        my $percent_gc = sprintf("%.2f", ($gc/$len) * 100);
        my $percent_at = sprintf("%.2f", ($at/$len) * 100);
        my $percent_ga = sprintf("%.2f", ($ga/$len) * 100);
        my $percent_ct = sprintf("%.2f", ($ct/$len) * 100);
        my $percent_gt = sprintf("%.2f", ($gt/$len) * 100);
        my $percent_ac = sprintf("%.2f", ($ac/$len) * 100);
        my $percent_nn = sprintf("%.2f", ($nn/$len) * 100);

        print TSV $basename."\t";
        print TSV $key."\t";
        print TSV $percent_gc."\t";
        print TSV $percent_at."\t";
        print TSV $percent_ga."\t";
        print TSV $percent_ct."\t";
        print TSV $percent_gt."\t";
        print TSV $percent_ac."\t";
        print TSV $percent_nn."\t";

        my $f_len = commify($len);
        my $f_nn = commify($nn);
        print TSV $f_len."\t";
        print TSV $f_nn."\n";

    }

    print TSV "\n";

    ## Summary

    my $tot_percent_gc = sprintf("%.2f", ($tot_gc/$tot_len) * 100);
    my $tot_percent_at = sprintf("%.2f", ($tot_at/$tot_len) * 100);
    my $tot_percent_ga = sprintf("%.2f", ($tot_ga/$tot_len) * 100);
    my $tot_percent_ct = sprintf("%.2f", ($tot_ct/$tot_len) * 100);
    my $tot_percent_gt = sprintf("%.2f", ($tot_gt/$tot_len) * 100);
    my $tot_percent_ac = sprintf("%.2f", ($tot_ac/$tot_len) * 100);
    my $tot_percent_nn = sprintf("%.2f", ($tot_nn/$tot_len) * 100);

    print SUM $basename."\t";
    print SUM $tot_percent_gc."\t";
    print SUM $tot_percent_at."\t";
    print SUM $tot_percent_ga."\t";
    print SUM $tot_percent_ct."\t";
    print SUM $tot_percent_gt."\t";
    print SUM $tot_percent_ac."\t";
    print SUM $tot_percent_nn."\t";

    my $f_tot_len = commify($tot_len);
    my $f_tot_nn = commify($tot_nn);
    print SUM $f_tot_len."\t";
    print SUM $f_tot_nn."\n";

}

### Subroutine(s)
sub commify {
	my $text = reverse $_[0];
	$text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
	return scalar reverse $text;
}