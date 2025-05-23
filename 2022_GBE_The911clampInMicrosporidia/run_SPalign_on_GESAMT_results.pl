#!/usr/bin/perl
## Pombert Lab 2022

my $name = "run_SPalign_on_GESAMT_results.pl";
my $version = "0.1a";
my $updated = "2022-03-09";

use strict;
use warnings;
use PerlIO::gzip;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	The purpose of this script is to to calculate the GDT score by running the SP-align structural alignment tool,
		an indicator of the strength a predicted protein model shares with the expiremental model it is attempting to represent.

USAGE

OPTION
-t (--tdfi)		3DFI output folder
-r (--rcsb)		Path to RCSB PDB structures
-o (--outdir)		Output directory (Default: SPalign_RESULTS)
EXIT

die("\n$usage\n") unless(@ARGV);

my ($script,$pipeline_dir) = fileparse($0);
$pipeline_dir =~ s/Homology_search/Misc_tools/;

my %Folds = ("ALPHAFOLD" => "ALPHAFOLD_3D_Parsed",
			 "RAPTORX" => "RAPTORX_3D/PDB"
);

my $tdfi;
my @rcsb;
my $outdir = "SPalign_RESULTS";

GetOptions(
	'-t|--tdif=s{1}' => \$tdfi,
	'-r|--rcsb=s@{1,}' => \@rcsb,
	'-o|--outdir=s' => \$outdir,
);

my $rcsb_temp_dir = "$outdir/tmp_rcsb";

unless(-d $outdir){
	make_path($rcsb_temp_dir,{mode=>0755});
}

my $datestring = localtime();
open LOG, ">", "$outdir/SPalign.log" or die("Unable to create file $outdir/SPalign.log: $!\n");
print LOG ("$0 \\\n-t $tdfi \\\n-r @rcsb \\\n-o $outdir\n\n");
print LOG ("Started on $datestring\n");

open IN, "<", "$tdfi/Homology/GESAMT/All_GESAMT_matches_per_protein.tsv" or die "Unable to open file $tdfi/Homology/GESAMT/All_GESAMT_matches_per_protein.tsv: $!\n";
my $total_alignments;
while(my $line = <IN>){
	unless(($line =~ /^###/) || ($line eq "")){
		$total_alignments++;
	}
}
close IN;

open IN, "<", "$tdfi/Homology/GESAMT/All_GESAMT_matches_per_protein.tsv" or die "Unable to open file $tdfi/Homology/GESAMT/All_GESAMT_matches_per_protein.tsv: $!\n";
open RAW, ">", "$outdir/SPalign_raw.tsv";
print RAW ("### Query\tPredictor\tRCSB Code\tChain\tRMSD\tnAlign\tGDT\n");

my $alignment_counter = 0;
while(my $line = <IN>){
	chomp($line);
	unless(($line =~ /^###/) || ($line eq "")){


		my ($query,$predictor,$rcsb_code,$chain,$qscore,$rmsd,$seq_id,$nAlign,$nRes,$rcsb_file,$annotation) = split("\t",$line);
		my ($rcsb_sub_folder) = lc($rcsb_code) =~ /\w(\w{2})\w/;

		## Checking in current and obsolete RCSB PDB folders
		my $rcsb_file_location;
		for my $rcsb_dir (@rcsb){
			if (-f "$rcsb_dir/$rcsb_sub_folder/$rcsb_file"){
				$rcsb_file_location = "$rcsb_dir/$rcsb_sub_folder/$rcsb_file";
				last;
			}
		}

		if($rcsb_file_location){
			system "clear";
			print("\n\tAligning $query to $rcsb_code\n");
			my $remaining = "." x (int((($total_alignments-$alignment_counter)/$total_alignments)*100));
			my $progress = "|" x (100-int((($total_alignments-$alignment_counter)/$total_alignments)*100));
			my $status = "[".$progress.$remaining."]";
			print "\n\t$status\t".($alignment_counter)."/$total_alignments\n\t";

			system "cp $rcsb_file_location $rcsb_temp_dir/$rcsb_file\n";

			if(-f "$pipeline_dir/split_PDB.pl"){
				system "$pipeline_dir/split_PDB.pl \\
						-p $rcsb_temp_dir/$rcsb_file \\
						-o $rcsb_temp_dir/tmp/ \\
						-e pdb
				";
			}
			else{
				print STDERR "[E] Cannot find $pipeline_dir/split_PDB.pl\n";
			}

			my $temp_file = "$rcsb_temp_dir/tmp/pdb".lc($rcsb_code)."/pdb".lc($rcsb_code)."_$chain.pdb";
			my $predicted_file_location = "$tdfi/Folding/".$Folds{$predictor}."/$query.pdb";

			my $SPalign_result = `SP-align $predicted_file_location $temp_file`;
			$alignment_counter++;
			
			system "rm -rf $rcsb_temp_dir/*";

			my @data = split("\n",$SPalign_result);

			my $RMSD;
			my $nAlign;
			my $GDT;
			my $TMscore;
			foreach my $line (@data){
				chomp($line);
				if($line =~ /GDT/){
					my ($RMSD_seg,$GDT_seg,$TMscore_seg) = split(";",$line);
					($RMSD,$nAlign) = $RMSD_seg =~ /^RMSD\/Nali= (\S+) \/ (\d+)/;
					($GDT) = $GDT_seg =~ /^GDT= (\S+)/;
					print RAW ("$query\t$predictor\t$rcsb_code\t$chain\t$RMSD\t$nAlign\t$GDT\n");
				}
			}
		}
		else{
			print STDERR "Unable to find $rcsb_file\n";
		}

	}
	elsif($line =~ /^### (\w+)/){
		print RAW ("### $1\n");
	}
	else{
		print RAW ("\n");
	}
}

print LOG ("$alignment_counter proteins aligned\n");
$datestring = localtime();
print LOG ("Completed on $datestring\n");

close IN;
close RAW;
close LOG;
