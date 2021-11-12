#!/usr/bin/perl
my $name = "add_b_values.pl";
my $version = "0.2.8";
my $updated = "2021-10-18";

use strict; use warnings; use File::Basename; use Getopt::Long qw(GetOptions); use Cwd qw(abs_path);

my $usage = << "EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Adds plDDT values to the B-factor column of alphafold predictions that were performed
		before this feature was implemented. Also extracts these scores as tab-delimited files from the
		corresponding .pkl files.  
EXAMPLE		${name} \\
		  -d ./Alphafold \\
		  -c AF_colored_PDB \\
		  -e AF_plDDT_scores
OPTIONS
-d (--dir)	Input directory
-c (--col)	Add confidence scores to AlphaFold .pdb files from .pkl file
-e (--ext)	Extract by-residue confidence scores from .pkl to indepent file
EXIT
die "\n$usage\n" unless @ARGV;

my $in_dir;
my $ext_dir = "AlphaFold_Extracted_Confidence_Scores";
my $col_dir = "AlphaFold_PDB_wConfidence_Scores";
GetOptions(
	"dir|d=s" => \$in_dir,
	"col|c=s" => \$col_dir,
	"ext|e=s" => \$ext_dir,
);

die "\n[E]  No process has been selected\n" unless ($col_dir||$ext_dir);

unless (-d $col_dir) { 
	mkdir ($col_dir, 0755) or die "Can't create $col_dir: $!\n";
}
unless (-d $ext_dir) { 
	mkdir ($ext_dir, 0755) or die "Can't create $ext_dir: $!\n";
}

## Converting input dir to absolute path
$in_dir = abs_path($in_dir);

## For each input file, check to see if it is a pickle file
opendir(EXTDIR,$in_dir) or die "Can't open $in_dir: $!\n";
foreach my $obj (readdir(EXTDIR)){

		if ($obj =~ /alphafold2.log/){ next; }
		elsif ($obj =~/^\./){ next; }
		elsif ($obj =~ /^msas/){ next; }
		elsif ($obj =~ /\.txt$/){ next; }

		## Working on each AlphaFold subfolder
		opendir(INTDIR,"$in_dir/$obj") or die "Can't open $in_dir/$obj: $!\n";;

		## Obtain alphafold ranking of models
		my %ranking;
		if (-f "$in_dir/$obj/ranking_debug.json"){
			open RANK,"<","$in_dir/$obj/ranking_debug.json" or die "Can't open $in_dir/$obj/ranking_debug.json: $!\n";
			my $record = 0;
			## Rank counter
			my $counter = 0;
			print($obj."\n");
			while (my $line = <RANK>){
				## Don't start recording model numbers until we have reached the ranking
				## section of the .json file
				if ($record == 1){
					if ($line =~ /"(model_\d+)"/){
						$ranking{$1} = $counter;
						$counter++;
					}
				}
				else{
					if ($line =~ /"order"/){
						$record = 1;
					}
				}
				
			}
		}
		## Prevent recursing into the output directory
		elsif ($obj =~ /$col_dir/){
			next;
		}
		else{
			print("\n\t[E]\t$in_dir/$obj does not contain a rank file. It may either be missing or this directory is not an AlphaFold output directory!\n\n");
			next;
		}
		foreach my $in_file (readdir(INTDIR)){
			if ($in_file =~ /.pkl$/){
				unless ($in_file =~ /features.pkl/){
					my $script_loc = abs_path($0);
					my ($script,$work_dir) = fileparse($script_loc);
					system("$work_dir/extract_b_values.py $in_dir/$obj/$in_file");
					open IN, "<", "$in_dir/$obj/$in_file.b_values" or die "Can't open $in_dir/$obj/$in_file.b_values: $!\n";
					my @b_values;
					while (my $b_val = <IN>){
						push(@b_values,$b_val);
					}
					my ($model_name) = $in_file =~ /result_(\w+).pkl/;
					close IN;
					my $model = "relaxed";
					open COL, ">", "$col_dir/${obj}-${model}_${model_name}_w_bvalue.pdb";
					open EXT, ">", "$ext_dir/${obj}-${model}_${model_name}_confidence_scores.tsv";
					open IN, "<", "$in_dir/$obj/${model}_${model_name}.pdb";
					## b-factor counter (matches atom to proper b-value)
					my $counter = -1;
					my $prev = 0;
					TEST: while (my $line = <IN>){
						chomp($line);
						## Check to see if the file contains regex-breaking "MODEL" line
						if ($line =~ /^MODEL/){
							next;
						}
						## Check to see if we have reached the end of the file
						if ($line =~ /^END/){
							last TEST;
						}
						my @data_set = split(/\s+/,$line);
						my $boundary_counter = 0;
						foreach my $data (@data_set){
							## PDB-format is special, so sometimes two values can get glued
							## together so check to see if these values need to be split
							if ($data =~ /(\D)(\d{4,})/){
								splice(@data_set,4,1,$1,$2);
							}
							## Three columns combined together
							if ($data =~ /(\-\S+)(\-\S+)(\-\S+)/){
								splice(@data_set,6,1,$1,$2,$3);
							}
							elsif ($data =~ /(\S+)(\-\S+)(\-\S+)/){
								splice(@data_set,6,1,$1,$2,$3);
							}
							## Two columns combined together
							elsif ($data =~ /(\-\S+)(\-\S+)/){
								if ($data_set[$boundary_counter+1] eq '1.00'){ ## right bounded
									splice(@data_set,7,1,$1,$2);
								}
								else { ## left bounded
									splice(@data_set,6,1,$1,$2);
								}
							}
							elsif ($data =~ /(\S+)(\-\S+)/){
								if ($data_set[$boundary_counter+1] eq '1.00'){ ## right bounded
									splice(@data_set,7,1,$1,$2);
								}
								else { ## left bounded
									splice(@data_set,6,1,$1,$2);
								}
							}

							$boundary_counter++;
							
						}
						## Print all the values in proper, PDB approved, format
						my $out_line = sprintf("%-6s",$data_set[0]);
						$out_line .= sprintf("%5d",$data_set[1]);
						$out_line .= " ";
						$out_line .= sprintf("%-4s",$data_set[2]);
						$out_line .= sprintf("%1s",'');
						$out_line .= sprintf("%3s",$data_set[3]);
						$out_line .= " ";
						$out_line .= sprintf("%1s",$data_set[4]);
						## If we have reached to terminate line, print only the existing values,
						## if not, print all values
						unless ($line =~ /TER/){
							unless($prev == $data_set[5]){
								$prev = $data_set[5];
								print EXT $data_set[5]."\t".$data_set[3]."\t".$b_values[$counter];
								$counter++;
							}
							$out_line .= sprintf("%4d",$data_set[5]);
							$out_line .= sprintf("%1s",'');
							$out_line .= "   ";
							$out_line .= sprintf("%8.3f",$data_set[6]);
							$out_line .= sprintf("%8.3f",$data_set[7]);
							$out_line .= sprintf("%8.3f",$data_set[8]);
							$out_line .= sprintf("%6.2f",$data_set[9]);
							$out_line .= sprintf("%6.2f",$b_values[$counter]);
							$out_line .= "          ";
							$out_line .= sprintf("%2s",$data_set[11]);
							$out_line .= sprintf("%2s",'');
						}
						## Print the line out to the new PDB file
						print COL $out_line."\n";
					}
					## Get the current model name so we can reference the proper rank from the
					## ranking_debug.json
					my $current_model = $model_name;
					## Copy the relaxed_w_bvalue pdb as ranked pdb in proper placement
					system ("cp $col_dir/${obj}-relaxed_${current_model}_w_bvalue.pdb $col_dir/${obj}-r".$ranking{$current_model}.".pdb");
					print ("  Fixed ${obj}-r$ranking{$current_model}\n");
					system ("rm -rf $col_dir/${obj}-relaxed_${current_model}_w_bvalue.pdb");
				}
			}
		}

}
