#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'run_megadock.pl';
my $version = '0.2';
my $updated = '2022-06-01';

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
SYNOPSIS	Runs megadock from its compiled binary or docker image

EXAMPLE		${name} \\
		  -o output_dir \\
		  -r receptors/*.pdb \\
		  -l ligands/*.pdb \\
		  -pr 3 \\
		  -po 10000

GENERAL OPTIONS:
-o (--outdir)		Output directory [Default: MDOCK]
-r (--receptors)	PDB files of proteins to be used as receptors
-l (--ligands)		PDB files of proteins to be used as ligands
--table		Tab-delimited list of receptors and ligands to dock

# Docking options
-pr (--pred_rot)	Number of predictions per each rotation [Default: 3]
-po (--pred_out)	Number of output predictions [Default: 10000]
-d (--detail)		Output docking detail files (.tsv and .detail)
--runmode		Runmode: [Default: megadock-gpu]
			  megadock-gpu		# GPU, single node
			  megadock-gpu-dp	# GPU, multiple nodes
			  megadock		# CPU, single node
			  megadock-dp		# CPU, multiple nodes
--nodes			Number of MPI nodes to use [Default: 10]
--host			Host file for MPI nodes

# Docker options
--docker		Uses the Megadock docker image instead of compiled version
--image			Docker image name [Default: akiyamalab/megadock:gpu]
USAGE
die "\n$usage\n" unless @ARGV;

my $outdir = 'MDOCK';
my @receptors;
my @ligands;
my $table;

my $pred_rot = 3;
my $pred_out = 10000;
my $detail;
my $runmode = 'megadock-gpu';
my $mpinodes = 10;
my $host;

my $docker;
my $image = 'akiyamalab/megadock:gpu';

GetOptions(
	'o|outdir=s' => \$outdir,
	'r|receptors=s@{0,}' => \@receptors,
	'l|ligands=s@{0,}' => \@ligands,
	'table=s' => \$table,
	'pr|pred_rot=i' => \$pred_rot,
	'po|pred_out=i' => \$pred_out,
	'd|detail' => \$detail,
	'runmode=s' => \$runmode,
	'nodes=i' => \$mpinodes,
	'host=s' => \$host,
	'docker' => \$docker,
	'image=s' => \$image,
);

##### Check for ligands and receptors
unless ($table){
	unless (@receptors){
		print "Please provide the receptors/ligands to dock with the --table or the -r and -l command line options.\n";
		exit;
	}
}


##### Grabbing user and group ids from the shell
# will use to run as user rather than root
my $uid = `id -u`;
my $gid = `id -g`;
chomp $uid;
chomp $gid;

##### Check if output directory / subdirs can be created
$outdir =~ s/\/$//;
unless (-d $outdir) {
	make_path($outdir, {mode => 0755 }) or die "Can't create $outdir: $!\n";
}
$outdir = abs_path($outdir);

# Create subdirs to store data
my $REC = "$outdir/RECEPTORS";		## Store proteins to be used as receptors
my $LIG = "$outdir/LIGANDS";		## Store proteins to be used as ligands
my $DOC = "$outdir/DOCKING";		## Store output of docking, create a subdir per protein receptor
for my $subdir ($REC,$LIG,$DOC){
	unless (-d $subdir) {
		make_path($subdir, {mode => 0755}) or die "Can't create $subdir: $!\n";
	}
}

## Setting data path for Megadock4 docker image
my $data_path = $outdir;

## Creating list of docking computations to perform
my %todock;
if ($table){
	open TABLE, "<", "$table" or die "Can't open $table: $!\n";
	while (my $line = <TABLE>){
		chomp $line;
		my @data = split ("\t", $line);
		my $rec = $data[0];
		my $lig = $data[1];
		push (@{$todock{$rec}}, $lig);
	}
}
if (@receptors){
	foreach my $rec (@receptors){
		foreach my $lig (@ligands){
			push (@{$todock{$rec}}, $lig);
		}
	}
}

my @recs = keys %todock;

## Iterating through all receptors

print "\n### Running Megadock 4\n";

foreach my $receptor (@recs){

	my ($rec_name) = fileparse($receptor);
	my ($rec_prefix) = $rec_name =~ /(\S+)\.pdb$/;
	my $filecount = 0;

	## Copying receptor to output directory
	unless (-e "$REC/$rec_name"){
		system "cp $receptor $REC/";
	}

	## Creating subdir to store docking results per protein receptor
	my $rec_dir = "$outdir/DOCKING/$rec_prefix";
	unless (-e "$rec_dir"){
		make_path( $rec_dir, { mode => 0755 } )  or die "Can't create $rec_dir: $!\n";
	}

	## Iterating through all ligands
	my @ligs = @{$todock{$receptor}};
	foreach my $ligand (@ligs){

		my ($lig_name) = fileparse($ligand);
		my ($lig_prefix) = $lig_name =~ /(\S+)\.pdb$/;
		$filecount++;

		## Copying ligand to output directory
		unless (-e "$LIG/$lig_name"){
			system "cp $ligand $LIG/";
		}

		## Setting filenames
		my $recfile = "$outdir/RECEPTORS/$rec_name";
		my $ligfile = "$outdir/LIGANDS/$lig_name";
		my $outfile = "$outdir/DOCKING/$rec_prefix/$rec_prefix.vs.$lig_prefix.out";

		## Adjusting filenames to account if --docker is requested
		if ($docker){
			$recfile = "data/RECEPTORS/$rec_name";
			$ligfile = "data/LIGANDS/$lig_name";
			$outfile = "data/DOCKING/$rec_prefix/$rec_prefix.vs.$lig_prefix.out";
		}

		my $realfile = "$outdir/DOCKING/$rec_prefix/$rec_prefix.vs.$lig_prefix.out";

		## Running docking with Megadock 4
		if (-e $realfile){
			my $modulo = $filecount % 1000;
			if ($modulo == 0){
				print "  $filecount: $realfile already exists. Skipping ...\n";
			}
		}
		else {

			## Check for MPI nodes
			my $mpirun_flag;
			if (($runmode eq 'megadock-gpu-dp') or ($runmode eq 'megadock-dp')){
				$mpirun_flag = 1;
			}

			# Check for MPI host file
			my $mpi_host_flag = '';
			if ($host){
				$mpi_host_flag = "-hostfile $host";
			}

			## Check for detail flag
			my $detail_flag = '';
			if ($detail){
				$detail_flag = '-O';
			}

			## Use docker image if --docker flag is set
			if ($docker){
				system ("docker \\
					run \\
					-it \\
					--user $uid:$gid \\
					--runtime=nvidia \\
					--privileged \\
					-v $data_path:/opt/MEGADOCK/data \\
					$image \\
					$runmode \\
					-R $recfile \\
					-L $ligfile \\
					-t $pred_rot \\
					-N $pred_out \\
					$detail_flag \\
					-o $outfile
				") == 0 or checksig();
			}

			## Use compiled binaries

			# mpirun
			elsif ($mpirun_flag){

				# Creating table file for MPI run
				my $table = "$outdir/DOCKING/$rec_prefix/$rec_prefix.vs.$lig_prefix.table";
				open TABLE, ">", $table or die "Can't create $table: $!\n";
				print TABLE "TITLE=$rec_prefix.vs.$lig_prefix\n";
				# $2 is not recognized better to enter the $ligfile directly
				print TABLE "PARAM=-R $recfile -L $ligfile $detail_flag -o $outfile\n";
				# Medadock expect a 3rd line; A commented line is fine
				print TABLE "$recfile $ligfile";

				system ("mpirun -n $mpinodes \\
					$mpi_host_flag \\
					$runmode \\
					-tb $table;
				") == 0 or checksig();
				## Somehow returns exit code 1 even when works...
				## MPI_ABORT is invoked, not invoking MPI_Finalize() properly?
			}

			# single node
			else {
				system ("$runmode \\
					-R $recfile \\
					-L $ligfile \\
					-t $pred_rot \\
					-N $pred_out \\
					$detail_flag \\
					-o $outfile
				") == 0 or checksig();
			}
		}

	}

}

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
