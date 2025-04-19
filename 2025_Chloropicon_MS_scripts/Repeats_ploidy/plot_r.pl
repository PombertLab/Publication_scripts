#!/usr/bin/env perl
## Pombert Lab, Illinois Tech, 2023
my $name = 'plot_r.pl';
my $version = '0.1';
my $updated = '2023-05-08';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd qw(abs_path);

my $usage =<<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Plots histograms with R

COMMAND		${name} \\
		  -t *.tsv \\
		  -o outdir \\
		  -ymax 0.12

OPTIONS:
-t (--tsv)	TSV files to plots
-o (--out)	Output directory
-xmin		Minimum x value [Default: 0]
-xmax		Maximum x value [Default: 100]
-ymin		Minimum y value [Default: 0]
-ymax		Maximum y value [Default: 0.12]
OPTIONS
die "\n$usage\n" unless @ARGV;

my @tsv;
my $outdir;
my $xmin = 0;
my $xmax = 100;
my $ymin = 0;
my $ymax = 0.12;
GetOptions(
	't|tsv=s@{1,}' => \@tsv,
	'o|out=s'	=> \$outdir,
	'xmin=s' => \$xmin,
	'xmax=s' => \$xmax,
	'ymin=s' => \$ymin,
	'ymax=s' => \$ymax
);

## Outdir
unless (-d $outdir){
	mkdir ($outdir, 0755) or die $!;
}

## Rscript
my $rscript = abs_path($outdir.'/'.'plots.R');

open R, '>', $rscript or die $!;
print R '#!/usr/bin/env Rscript'."\n\n";

while (my $file = shift @tsv){
	rprint($file,$outdir,\*R);
}

close R;

system "chmod +x $rscript";
system "$rscript";

### Subroutine
sub rprint {

	my $file = abs_path($_[0]);
	my $dir = abs_path($_[1]);
	my $fh = $_[2];
	my ($basename) = fileparse($file);
	my ($prefix) = $basename =~ /^(\S+)\.\w+/;

	print $fh 'matrix <- as.matrix(read.table("'.$file.'", header=FALSE, sep="\t"))'."\n";
	print $fh 'pdf("'.$outdir.'/'.$prefix.'.pdf'.'")'."\n";

	print $fh 'library(MASS)'."\n";
	print $fh 'colnames(matrix)[1] <- "allelic frequency"'."\n";
	print $fh 'fit <- fitdistr(matrix, "normal")'."\n";
	print $fh 'class(fit)'."\n";
	print $fh 'para <- fit$estimate'."\n";

	print $fh 'hist(matrix, breaks = 50, prob = TRUE, xlim=c('.$xmin.','.$xmax.'), ylim=c('.$ymin.','.$ymax.'), xlab = "allelic frequencies")'."\n";
	print $fh 'curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)'."\n";

	print $fh 'while (!is.null(dev.list()))  dev.off()'."\n\n";

}