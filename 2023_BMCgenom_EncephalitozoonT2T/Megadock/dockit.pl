#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'dockit.pl';
my $version = '0.2a';
my $updated = '2022-09-29';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);

my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs protein-protein docking with Megadock 4, get PPIscores and
		generates PDB files from the best matches

EXAMPLE		${name} \\
		  -o MDOCK \\
		  -r *.pdb \\
		  -l *.pdb \\
		  -p Eint_products.tsv 

OPTIONS:
### General
-o (--outdir)		Output directory [Default: MDOCK]
-r (--receptors)	PDB files of proteins to be used as receptors
-l (--ligands)		PDB files of proteins to be used as ligands
-p (--products)	Tab-delimited list of ligands and their products
--table			Tab-delimited list of receptors and ligands to dock

### Docking
-pr (--pred_rot)	Number of predictions per each rotation [Default: 3]
-po (--pred_out)	Number of output predictions [Default: 10000]
-d (--detail)		Output docking detail files (.tsv and .detail)
--nodes			Number of MPI nodes to use [Default: 10]
--host			Host file for MPI nodes
--runmode		Runmode: [Default: megadock-gpu]
			  megadock-gpu		# GPU, single node
			  megadock-gpu-dp	# GPU, multiple nodes
			  megadock		# CPU, single node
			  megadock-dp		# CPU, multiple nodes
-ndo (--nodock)	Skip docking and go straight to PPI scoring

### PPI/PDB
-ppi (--ppicut)		Minimum PPIscore to keep [Default: 5]
-t (--topmatch)		Maximum number of matches to keep [Default: 100]
-n (--number)		Number of PDBs to generate from top decoys [Default: 5]
-cat (--concat)		Concatenate all decoys into one PDB file # Instead of one decoy per PDB

### Docker
--docker		Uses the Megadock docker image instead of compiled version
--image			Docker image name [Default: akiyamalab/megadock:gpu]
USAGE
die "\n$usage\n" unless @ARGV;

## General
my $outdir = 'MDOCK';
my @receptors;
my @ligands;
my $products;
my $table;

## Docking
my $pred_rot = 3;
my $pred_out = 10000;
my $detail;
my $runmode = 'megadock-gpu';
my $mpinodes = 10;
my $host;
my $nodock;

## PPIscores
my $ppicutoff = 5;
my $tsvfile = 'ppiscores.tsv';
my $top_matches = 100;
my $num_decoys = 5;
my $concatenate;

## Docker
my $docker;
my $image = 'akiyamalab/megadock:gpu';

GetOptions(
	# General
	'o|outdir=s' => \$outdir,
	'r|receptors=s@{1,}' => \@receptors,
	'l|ligands=s@{1,}' => \@ligands,
	'p|products=s' => \$products,
	'table=s' => \$table,
	# Docking
	'pr|pred_rot=i' => \$pred_rot,
	'po|pred_out=i' => \$pred_out,
	'd|detail' => \$detail,
	'runmode=s' => \$runmode,
	'nodes=i' => \$mpinodes,
	'host=s' => \$host,
	'ndo|nodock' => \$nodock,
	# PPIscores
	'ppi|ppicut=i' => \$ppicutoff,
	't|topmatch=i' => \$top_matches,
	'n|number=i' => \$num_decoys,
	'cat|concat' => \$concatenate,
	# Docker
	'docker' => \$docker,
	'image=s' => \$image,
);

## Check output directory
$outdir =~ s/\/$//;
unless (-d $outdir) {
	make_path($outdir, {mode => 0755 }) or die "Can't create $outdir: $!\n";
}

########################################################################################
## Run Megadock 4
MEGADOCK:

if ($nodock){
	print "\nSkipping docking as requested. Moving to PPI scoring...\n\n";
	goto PPISCORES;
}

# Check for detail flag
my $detail_flag = '';
if ($detail){
	$detail_flag = '-O';
}

# Check for MPI host file
my $mpi_host_flag = '';
if ($host){
	$mpi_host_flag = "--host $host";
}

# Check for Docker flag
my $docker_flag = '';
my $docker_img_flag = '';
if ($docker){
	$docker_flag = "--docker";
	$docker_img_flag = "--image $image";
}

### Receptors and ligands
my $todock_list = $outdir.'/'.'todock.tsv';

# Use table if provided
if ($table){
	system "cp $table $outdir/todock.tsv";
}
# Otherwise create one from -r + -l
else {
	open TODO, ">", $todock_list or die "Can't create $todock_list: $!\n";
	foreach my $currec (@receptors){
		foreach my $curlig (@ligands){
			print TODO "$currec\t$curlig\n";
		}
	}
	close TODO;
}

# Docking
system ("run_megadock.pl \\
	--table $todock_list \\
	--outdir $outdir \\
	--pred_rot $pred_rot \\
	--pred_out $pred_out \\
	--runmode $runmode \\
	$detail_flag \\
	$mpi_host_flag \\
	$docker_flag \\
	$docker_img_flag
") == 0 or checksig(); 

##### PPIscores
PPISCORES:

my $prod_flag = '';
if ($products){
	$prod_flag = "--products $products";
}

system ("get_ppiscores.pl \\
	--mdock $outdir \\
	--outdir $outdir/PPIs \\
	--npred $pred_out \\
	--tsv $tsvfile \\
	--cutoff $ppicutoff \\
	$prod_flag \\
	--verbose
	") == 0 or checksig();

##### PDBs
PDBS:

# Check for concatenate flag
my $concat_flag = '';
if ($concatenate){
	$concat_flag = '--concat';
}

system ("get_decoy_pdbs.pl \\
	--datadir $outdir \\
	--list $outdir/PPIs/$tsvfile \\
	--outdir $outdir/PDBs \\
	--minppi $ppicutoff \\
	--number $num_decoys \\
	$concat_flag
	") == 0 or checksig();

########################################################################################
## Exit

print "\nAll done. Exiting...\n\n";
exit;

########################################################################################
## Subroutine(s)
sub checksig {

	my $exit_code = $?;
	my $modulo = $exit_code % 255;

	print "\nExit code = $exit_code; modulo = $modulo \n";

	if ($modulo == 2) {
		print "\nSIGINT detected: Ctrl+C => exiting...\n\n";
		exit(2);
	}
	elsif ($modulo == 131) {
		print "\nSIGTERM detected: Ctrl+\\ => exiting...\n\n";
		exit(131);
	}

}
