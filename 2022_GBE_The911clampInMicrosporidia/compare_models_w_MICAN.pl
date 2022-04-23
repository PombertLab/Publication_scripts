#!/usr/bin/perl
## Pombert Lab 2022
my $name = 'compare_models_w_MICAN.pl';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

die("\n$name -t 3DFI_Results -o AF_compared_w_MICAN.tsv\n") unless(@ARGV);

my $tdfi_dir;
my $outfile = "AF_compared_w_MICAN.tsv";

GetOptions(
	't|tdfi=s' => \$tdfi_dir,
	'o|outfile=s' => \$outfile,
);

my $alpha_dir = "$tdfi_dir/Folding/ALPHAFOLD_3D/";

opendir(ALPHA,$alpha_dir) or die("Unable to open $alpha_dir: $!\n");

my %results;

foreach my $locus (readdir(ALPHA)){
	# print("$locus\n");
	my @ranked_files;
	my $min_pLDDT;
	my $max_pLDDT;
	unless(($locus eq ".") or ($locus eq "..")){
		if(-d "$alpha_dir/$locus"){
				opendir(LOCUS,"$alpha_dir/$locus");
				foreach my $ranked (readdir(LOCUS)){
					my $pdb = "$alpha_dir/$locus/$ranked";
					unless(-d $pdb){
						if($ranked =~ /ranked_\d\.pdb/){
							push(@ranked_files,$pdb);
						}
						if($ranked eq "ranking_debug.json"){
							open RANK, "<", $pdb;
							while(my $line = <RANK>){
								chomp($line);
								if ($line =~ /\"model_\d\": (.*)/){
									my $pLDDT = $1;
									$pLDDT =~ s/,//;
									$pLDDT = sprintf("%.3f",$pLDDT);
									if(!defined $min_pLDDT){
										$min_pLDDT = $pLDDT;
										$max_pLDDT = $pLDDT;
									}
									elsif($min_pLDDT > $pLDDT){
										$min_pLDDT = $pLDDT;
									}
									elsif($max_pLDDT < $pLDDT){
										$max_pLDDT = $pLDDT;
									}

								}
							}
							close RANK;
						}
					}
				}
			}

		print("Running MICAN now on $locus\n");

		## Mican
		my @mican_results;
		
		my $min;
		my $max;
		my $sum;
		my $denominator;

		my %permutations;

		foreach my $file1 (@ranked_files){
			foreach my $file2 (@ranked_files){
				unless($file1 eq $file2){
					my ($r1) = $file1 =~ /(\d)\.pdb$/;
					my ($r2) = $file2 =~ /(\d)\.pdb$/;
					my $mican_result = `mican -s $file1 $file2 -n 1`;
					my @data = split("\n",$mican_result);
					foreach my $tmp (@data){
						chomp $tmp;
						if ($tmp =~ /^\s+(1.*)/){
							my ($rank,$sTMscore,$TMscore,$Dali_Z,$SPscore,$Length,$RMSD,$Seq_Id) = split(/\s+/,$1);
							if (!defined $min){
								$min = $TMscore;
								$max = $TMscore;
							}
							elsif ($TMscore < $min){
								$min = $TMscore;
							}
							elsif ($TMscore > $max){
								$max = $TMscore;
								
							}
							$sum += $TMscore;
							$denominator++;
						}
					}
				}
			}
		}

		if ($denominator){
			my $result_line = $min."\t".$max;
			my $average = $sum/$denominator;
			$average = sprintf("%.3f", $average);
			$result_line .= "\t".$average;
			$result_line .= "\t".$min_pLDDT."\t".$max_pLDDT;
			$results{$locus} = $result_line;
		}
		else {
			$results{$locus} = ('N/A'."\t"x5)."N/A";
		}

	}
}

## Printing out model comparison statistics in locus order
open OUT, ">", "$outfile" or die("Unable to write to $outfile: $!\n");
print OUT "## TM- and pLDDT scores per locus\n";
print OUT "## Locus\tMin TMs\tMax TMs\tAvg TMs\tMin pLDDT\tMax pLDDT\n";
foreach my $key (sort(keys(%results))){
	print OUT "$key\t$results{$key}\n";
}