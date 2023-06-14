#!/usr/bin/perl
## Pombert Lab, Illinois Tech 2022

my $name = 'b2links.pl';
my $version = '0.1';
my $updated = '2022-05-25';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Parses the output of blastn searches (outfmt 6) of a genome against
		itself to create links for Circos

COMMAND		${name} \\
		  -b *.blastn.6 \\
		  -o Parsed/ \\
		  -m 5000 \\
		  -c

OPTIONS:
-b (--blast)		Blastn outfmt 6 file(s) to parse
-o (--outdir)		Output directory [Default: ./]
-mb (--minbit)		Minimum bitscore to keep [Default: 5000]
-ml (--minlen)		Minimum length to keep [Default: 1000]
-c (--circos)		Create a Circos links file
USAGE
die "\n$usage\n" unless @ARGV;

my @blast;
my $outdir = './';
my $minbit = 5000;
my $minlen = 1000;
my $circos;
GetOptions(
	'b|blast=s@{1,}' => \@blast,
	'o|outdir' => \$outdir,
	'mb|minbit=i' => \$minbit,
	'ml|minlen=i' => \$minlen,
	'c|circos'	=> \$circos
);

### Check if output directory / subdirs can be created
$outdir =~ s/\/$//;
unless (-d $outdir) {
	make_path($outdir,{mode => 0755})  or die "Can't create $outdir: $!\n";
}

### Iterating through BLAST files

while (my $blast = shift@blast){

	my $filename = fileparse($blast);
	my ($prefix) = $filename =~ /^(\w+)/;
	my $outfile = $outdir.'/'.$prefix.'.parsed.tsv';
	my $outlinks = $outdir.'/'.$prefix.'.links';

	open BLAST, "<", $blast or die "Can't open $blast: $!\n";
	open OUT, ">", $outfile or die "Can't create $outfile\n";
	
	if ($circos){
		open LINKS, ">", $outlinks or die "Can't create $outlinks\n";
		print LINKS "#locus1 start end locus2 start end\n";
	}

	my %queries;

	while (my $line = <BLAST>){

		chomp $line;
		my @data = split ("\t", $line);
		
		my $qseqid = $data[0];
		my $subject = $data[1];
		my $pident = $data[2];
		my $length = $data[3];
		my $mismatch = $data[4];
		my $gapopen = $data[5];
		my $qstart = $data[6];
		my $qend = $data[7];
		my $sstart = $data[8];
		my $send = $data[9];
		my $evalue = $data[10];
		my $bitscore = $data[11];

		## Adjustments for Circos
		my $s1 = $qstart - 1;
		my $s2 = $sstart - 1;
		my $e1 = $qend - 1;
		my $e2 = $send - 1;

		## check for self hits
		if ($qseqid eq $subject){
			if ($pident >= 100){
				if (($qstart == $sstart) and ($qend == $send)){
					## Self hit on full sequence
					next;
				}
			}
			else {
				if ($qstart > $sstart){
					## inverted hit, keepf the forward hit only
					next;
				}
				else {
					if ($bitscore >= $minbit){
						if ($length >= $minlen){
							print OUT "$line\n";
							if ($circos){
								print LINKS "$qseqid $s1 $e1 $subject $s2 $e2\n";
							}
						}
					}
				}
			}
		}

		## not a hit to self
		else {
			if ($bitscore >= $minbit){
				if ($length >= $minlen){
					print OUT "$line\n";
					if ($circos){
						print LINKS "$qseqid $s1 $e1 $subject $s2 $e2\n";
					}
				}
			}
		}
	
	}

}