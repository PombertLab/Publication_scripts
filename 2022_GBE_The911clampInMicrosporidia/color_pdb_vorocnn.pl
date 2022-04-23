#!/usr/bin/perl
my $name = 'color_pdb_vorocnn.pl';
my $version = '0.1a';
my $updated = '2021-09-01';

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPIS		Adds VoroCNN residue-level quality assessment to b-factor columns
NOTE		VoroCNN - https://team.inria.fr/nano-d/software/vorocnn/
EXAMPLE		${name} \\
		  -p PDB/ \\
		  -v vorocnn_output/ \\
		  -o PDB_voro_colored \\
		  -n
OPTIONS:
-p (--pdb)	PDB input folder
-v (--voro)	Voro quality scores folder
-o (--out)	Output folder [Default: PDB_voro_colored]
-n (--norm)	Normalize VoroCNN values so that max value = 100%
OPTIONS
die "\n$usage\n" unless @ARGV;

my $voro_dir;
my $pdb_dir;
my $outdir = 'PDB_voro_colored';
my $normalize;
GetOptions(
	'v|voro=s' => \$voro_dir,
	'p|pdb=s' => \$pdb_dir,
	'o|out=s' => \$outdir,
	'n|norm' => \$normalize
);

## Checking output directory
unless (-d $outdir){ mkdir ($outdir, 0755) or die "Can't create $outdir: $!\n"; }

## Reading VoroCNN quality scores from input folder
opendir (VORODIR, $voro_dir) or die "Can't open VoroCNN input directory $voro_dir: $!\n";
my @voro_scores;
while (my $vorocnn = readdir(VORODIR)){
	if ($vorocnn =~ /.scores$/){ 
		push (@voro_scores, "$voro_dir/$vorocnn");
	}
}
@voro_scores = sort @voro_scores;
closedir VORODIR;

## Working on individual VoroCNN files and corresponding PDBs
while (my $voro_score = shift @voro_scores) {

	## Pass #1; checking for max value
	open VORO, "<", "$voro_score" or die "Can't open $voro_score: $!\n";
	my %voro_db;
	my $max = 0;
	while (my $line = <VORO>){
		if ($line =~ /^chain_id/){ next; }
		else {
			my @columns = split("\t", $line);
			my $score = $columns[3];
			if ($score > $max){
				$max = $score;
			}
		}
	}
	close VORO;

	## Pass #2; working on data
	open VORO, "<", "$voro_score" or die "Can't open $voro_score: $!\n";
	while (my $line = <VORO>){
		chomp $line;
		if ($line =~ /^chain_id/){ next; }
		else {
			my @columns = split("\t", $line);
			my $position = $columns[2];
			my $score = $columns[3];
			if ($normalize){ ## Normalizing so that max = 100%, easier to compare accross figures (same legend)
				my $norm_score = ($score/$max)*100;
				$norm_score = sprintf("%.2f", $norm_score);
				if ($norm_score eq '100.00'){ $norm_score = '100.0'; }
				$voro_db{$position} = $norm_score;
			}
			else {
				$score = sprintf("%.2f", $score);
				$voro_db{$position} = $score;
			}
		}
	}

	## Working on PDB file
	my $pdb = fileparse($voro_score);
	$pdb =~ s/\.scores$//;
	$pdb = "$pdb_dir/$pdb.pdb";
	open PDB, "<", "$pdb" or die "Can't open $pdb: $!\n";
	my ($filename) = fileparse($pdb);

	open OUT, ">", "$outdir/$filename" or die "Can't create $outdir/$filename: $!\n";
	while (my $line = <PDB>){
		chomp $line;
		unless ($line =~ /^ATOM/){ print OUT "$line\n"; }
		else {
			if ($line =~ /(^ATOM\s+\d+\s+\S+(?:\s+)?(\w+)\s+(\w+\s+)?(\d+)\s+.*1.00)\s*\S+(\s+\w+\s{0,})$/){
				my $left_b = $1;
				my $aa = $2;
				my $chain = $3;
				my $position = $4;
				my $right_b = $5;
				print OUT "$left_b"." $voro_db{$position}"."$right_b\n";
			}
		} 
	}

}