#!/usr/bin/perl
## Pombert Lab 2022

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = <<"EXIT";
$0 -a product_file.tsv -p plddt_file.tsv -g gesamt_file.tsv -m mican_file.tsv -v vorocnn_file.tsv
EXIT

die("\n$usage\n") unless(@ARGV);

my $ncbi_file;
my $msdb_file;
my $plddt_file;
my $gesamt_file;
my $mican_file;
my $vorocnn_file;
my $spalign_file;

GetOptions(
	'n|ncbi=s' => \$ncbi_file,
	'd|msdb=s' => \$msdb_file,
	'p|plddt=s' => \$plddt_file,
	'g|gesamt=s' => \$gesamt_file,
	'm|mican=s'=> \$mican_file,
	'v|vorocnn=s' => \$vorocnn_file,
	's|spalign=s' => \$spalign_file,
);

## The annotation file was created by parsing the gbff file
## The contents of the file is as follows:
## Locus tag \t Annotation
open ANNOT, "<", $ncbi_file or die("Can't open file $ncbi_file: $!\n");
my %ncbi_annotations;
while(my $line = <ANNOT>){
	chomp($line);
	unless(($line =~ /^###/)||($line eq "")){
		my ($locus,$annotation) = split("\t",$line);
		$locus =~ s/\D$//;
		$ncbi_annotations{$locus} = $annotation;
	}
}
close ANNOT;

open ANNOT, "<", $msdb_file or die("Can't open file $msdb_file: $!\n");
my %msdb_annotations;
while(my $line = <ANNOT>){
	chomp($line);
	unless(($line =~ /^###/)||($line eq "")){
		my ($locus,$annotation) = split("\t",$line);
		$locus =~ s/\D$//;
		$msdb_annotations{$locus} = $annotation;
	}
}
close ANNOT;

## The gesamt file was created by the parse_all_models_by_Q.pl
## The contents of the file is as follows:
## Locus tag (with model number) \t Predictor \t RCSB code \t Chain \t Q-score \t RMSD \t SeqId \t Nalign \t nRes \t RCSB file \t Description
open GESAMT, "<", $gesamt_file or die("Can't open file $gesamt_file: $!\n");
my %alignment_scores;
my $alphafold_counter = 1;
my $raptorx_counter = 1;
my $previous_locus = '';
while(my $line = <GESAMT>){
	chomp($line);
	unless(($line =~ /^###/)||($line eq '')){
		my ($locus_tag,$predictor,$rcsb_code,$rcsb_chain,$qscore,$rmsd,undef,$nalign,undef,undef,$description) = split("\t",$line);
		my ($locus) = $locus_tag =~ /(\w+)-m1/;
		unless($previous_locus eq $locus){
			$alphafold_counter = 1;
			$raptorx_counter = 1;
		}
		if(($predictor eq 'ALPHAFOLD')&&($alphafold_counter < 6)){
			push(@{$alignment_scores{$locus}{$predictor}{"qscore"}},[$qscore,$description,$rcsb_code,$rcsb_chain,$rmsd,$nalign]);
			$alphafold_counter++;
		}
		if(($predictor eq 'RAPTORX')&&($raptorx_counter < 6)){
			push(@{$alignment_scores{$locus}{$predictor}{"qscore"}},[$qscore,$description,$rcsb_code,$rcsb_chain,$rmsd,$nalign]);
			$raptorx_counter++;
		}
		$previous_locus = $locus;
	}
}

my %model_scores;

## The plddt file was created by parsing the json files
## The contents of the file is a follows:
## Locus tag \t Model \t plddt \t Ranked file \t Unranked file
open PLDDT, "<", $plddt_file or die("Can't open file $plddt_file: $!\n");
while (my $line = <PLDDT>){
	chomp($line);
	my ($locus,$pdb,$plddt_score,$ranked,$unranked) = split("\t",$line);
	# print("$locus\n");
	if($line =~ /\w+\-m1\.pdb/){
		# print("$line\n");
		$plddt_score = sprintf ("%.4f", $plddt_score);
		$model_scores{$locus}{"plddt"} = $plddt_score;
	}
}
close PLDDT;

## The mican file was created by run_MICAN_on_gesamt_matches.pl
## The content of the file is as follows:
## Locus tag (with model number) \t Predictor \t RCSB code \t Chain \t sTMscore \t TMscore \t DaliZ \t SPscore \t Length \t RMSD \t SeqId
open MICAN, "<", $mican_file or die("Can't open file $mican_file: $!\n");
while(my $line=<MICAN>){
	chomp($line);
	unless(($line =~ /^###/)||($line eq '')){
		my @data = split("\t",$line);
		my ($locus_tag,$predictor,$rcsb_code,$rcsb_chain,$stmscore,$tmscore,$daliz,$spscore,$length,$rmsd,$seqid) = split("\t",$line);
		my ($locus) = $locus_tag =~ /(\w+)-m\d/;
		$locus =~ s/\D$//;
		if($alignment_scores{$locus}{uc($predictor)}{"qscore"}){
			foreach my $alignment (@{$alignment_scores{$locus}{uc($predictor)}{"qscore"}}){
				# print(join("\t",@{$alignment})."\n");
				my ($qscore,$description,$g_rcsb_code,$g_rcsb_chain) = @{$alignment};
				if(($g_rcsb_code eq $rcsb_code)&&($g_rcsb_chain eq $rcsb_chain)){
					# print("Match found\n");
					push(@{$alignment_scores{$locus}{uc($predictor)}{"tmscore"}},[$tmscore,$rmsd,$length,$seqid]);
					# print(join("\t",@{$alignment_scores{$locus}{uc($predictor)}{"tmscore"}[0]})."\n");
				}
			}
		}
	}
}

## The mican file was created by run_MICAN_on_gesamt_matches.pl
## The content of the file is as follows:
## Locus tag (with model number) \t Predictor \t RCSB code \t Chain \t sTMscore \t TMscore \t DaliZ \t SPscore \t Length \t RMSD \t SeqId
open SPALIGN, "<", $spalign_file or die("Can't open file $spalign_file: $!\n");
while(my $line=<SPALIGN>){
	chomp($line);
	unless(($line =~ /^###/)||($line eq '')){
		my @data = split("\t",$line);
		my ($locus_tag,$predictor,$rcsb_code,$rcsb_chain,$rmsd,$nalign,$gdt) = split("\t",$line);
		my ($locus) = $locus_tag =~ /(\w+)-m\d/;
		$locus =~ s/\D$//;
		if($alignment_scores{$locus}{uc($predictor)}{"qscore"}){
			foreach my $alignment (@{$alignment_scores{$locus}{uc($predictor)}{"qscore"}}){
				# print(join("\t",@{$alignment})."\n");
				my ($qscore,$description,$g_rcsb_code,$g_rcsb_chain) = @{$alignment};
				if(($g_rcsb_code eq $rcsb_code)&&($g_rcsb_chain eq $rcsb_chain)){
					# print("Match found\n");
					push(@{$alignment_scores{$locus}{uc($predictor)}{"gdt"}},[$gdt,$rmsd,$nalign]);
					# print(join("\t",@{$alignment_scores{$locus}{uc($predictor)}{"tmscore"}[0]})."\n");
				}
			}
		}
	}
}

## The vorocnn file was created by vorocnn_average.pl
## The content of the file is as follows:
## Locus tag \t AlphaFold score \t RaptorX score
open VOROCNN, "<", $vorocnn_file or die("Can't open file $vorocnn_file: $!\n");
while(my $line = <VOROCNN>){
	chomp($line);
	unless(($line =~ /^Locus_tag/)||($line eq '')){
		# print("$line\n");
		$line =~ s/N\/A/---/g;
		my ($locus,$voro_alpha,$voro_raptor) = split("\t",$line);
		$model_scores{$locus}{"vorocnn"}{"ALPHAFOLD"} = $voro_alpha;
		$model_scores{$locus}{"vorocnn"}{"RAPTORX"} = $voro_raptor;
	}
}

my $line;
print("Locus\tMicrosporidiaDB\tNCBI\tpLDDT");
print(("\tVoroCNN\tRCSB Code\tRCSB Chain\tQ-score\tRMSD\t3D Alignment Length\tTMscore\tRMSD\t3D Alignment Length\tGDT\tRMSD\t3D Alignment Length\tDescription"x2)."\n");
ANNOT:foreach my $locus (sort(keys(%msdb_annotations))){
	for (my $i = 0; $i < 6; $i++){
		$line = '';
		my $skip = 0;
		if($i == 0){
			
			## Add annotation from MicrosporidiaDB
			$line .= $locus."\t".$msdb_annotations{$locus};

			if($ncbi_annotations{$locus}){
				$line .= "\t".$ncbi_annotations{$locus};
			}
			else{
				$line .= "\t"."---";
			}

			## pLDDT scores
			if($model_scores{$locus}{"plddt"}){
				$line .= "\t".$model_scores{$locus}{"plddt"};
			}
			else{
				$line .= "\t"."---";
			}

			## VoroCNN scores
			if($model_scores{$locus}{"vorocnn"}){
				my $voro_alpha = $model_scores{$locus}{"vorocnn"}{"ALPHAFOLD"};
				if(defined $voro_alpha){
					$line .= "\t".$voro_alpha;
				}
				else{
					$line .= "\t"."---";
				}
			}
			else{
				$line .= "\t"."---";
			}

			if($alignment_scores{$locus}{"ALPHAFOLD"}{"qscore"}){
				my ($qscore,$description,$rcsb_code,$rcsb_chain,$rmsd,$nalign) = @{$alignment_scores{$locus}{"ALPHAFOLD"}{"qscore"}[0]};
				$line .= "\t".$rcsb_code."\t".$rcsb_chain."\t".$qscore."\t".$rmsd."\t".$nalign;
				if($alignment_scores{$locus}{"ALPHAFOLD"}{"tmscore"}){
					my ($tmscore,$rmsd,$length,$seqid) = @{$alignment_scores{$locus}{"ALPHAFOLD"}{"tmscore"}[0]};
					# print("$locus\t$rcsb_code\t$rcsb_chain\t$tmscore\t$rmsd\n");
					$line .= "\t".$tmscore."\t".$rmsd."\t".$length;
				}
				else{
					$line .= "\t---"x3;
				}
				if($alignment_scores{$locus}{"ALPHAFOLD"}{"gdt"}){
					my ($gdt,$rmsd,$nalign) = @{$alignment_scores{$locus}{"ALPHAFOLD"}{"gdt"}[0]};
					# print("$locus\t$rcsb_code\t$rcsb_chain\t$tmscore\t$rmsd\n");
					$line .= "\t".$gdt."\t".$rmsd."\t".$nalign;
				}
				else{
					$line .= "\t---"x3;
				}
				$line .= "\t".$description;
			}
			else{
				$line .= "\t---"x12;
				$skip++;
			}

			## RaptorX segment of table
			if($model_scores{$locus}{"vorocnn"}){
				my $voro_raptor = $model_scores{$locus}{"vorocnn"}{"RAPTORX"};
				if(defined $voro_raptor){
					$line .= "\t".$voro_raptor;
				}
				else{
					$line .= "\t"."---";
				}
			}
			else{
				$line .= "\t"."---";
			}

			if($alignment_scores{$locus}{"RAPTORX"}{"qscore"}){
				my ($qscore,$description,$rcsb_code,$rcsb_chain,$rmsd,$nalign) = @{$alignment_scores{$locus}{"RAPTORX"}{"qscore"}[0]};
				$line .= "\t".$rcsb_code."\t".$rcsb_chain."\t".$qscore."\t".$rmsd."\t".$nalign;

				# @{$alignment_scores{$locus}{uc($predictor)}{"tmscore"}},[$tmscore,$rmsd])

				if($alignment_scores{$locus}{"RAPTORX"}{"tmscore"}){
					my ($tmscore,$rmsd,$length,$seqid) = @{$alignment_scores{$locus}{"RAPTORX"}{"tmscore"}[0]};
					# print("$locus\t$tmscore\t$rmsd\n");
					$line .= "\t".$tmscore."\t".$rmsd."\t".$length;
				}
				else{
					$line .= "\t---"x3;
				}
				if($alignment_scores{$locus}{"RAPTORX"}{"gdt"}){
					my ($gdt,$rmsd,$nalign) = @{$alignment_scores{$locus}{"RAPTORX"}{"gdt"}[0]};
					# print("$locus\t$rcsb_code\t$rcsb_chain\t$tmscore\t$rmsd\n");
					$line .= "\t".$gdt."\t".$rmsd."\t".$nalign;
				}
				else{
					$line .= "\t---"x3;
				}
				$line .= "\t".$description;
			}
			else{
				$line .= "\t---"x12;
				$skip++;
			}

			print("$line\n");
			if($skip >= 2){
				print("\n");
				next ANNOT;
			}

		}
		else{
			if($alignment_scores{$locus}{"ALPHAFOLD"}{"qscore"}[$i]){
				$line .= "\t"x4;
				my ($qscore,$description,$rcsb_code,$rcsb_chain,$rmsd,$nalign) = @{$alignment_scores{$locus}{"ALPHAFOLD"}{"qscore"}[$i]};
				$line .= "\t".$rcsb_code."\t".$rcsb_chain."\t".$qscore."\t".$rmsd."\t".$nalign;
				if($alignment_scores{$locus}{"ALPHAFOLD"}{"tmscore"}[$i]){
					my ($tmscore,$rmsd,$length,$seqid) = @{$alignment_scores{$locus}{"ALPHAFOLD"}{"tmscore"}[$i]};
					$line .= "\t".$tmscore."\t".$rmsd."\t".$length;
				}
				else{
					$line .= "\t"x6;
				}
				if($alignment_scores{$locus}{"ALPHAFOLD"}{"gdt"}[$i]){
					my ($gdt,$rmsd,$nalign) = @{$alignment_scores{$locus}{"ALPHAFOLD"}{"gdt"}[$i]};
					$line .= "\t".$gdt."\t".$rmsd."\t".$nalign;
				}
				else{
					$line .= "\t"x6;
				}
				$line .= "\t".$description;
			}
			else{
				$line .= "\t"x16;
				$skip++;
			}

			if($alignment_scores{$locus}{"RAPTORX"}{"qscore"}[$i]){
				$line .= "\t";
				my ($qscore,$description,$rcsb_code,$rcsb_chain,$rmsd,$nalign) = @{$alignment_scores{$locus}{"RAPTORX"}{"qscore"}[$i]};
				$line .= "\t".$rcsb_code."\t".$rcsb_chain."\t".$qscore."\t".$rmsd."\t".$nalign;
				if($alignment_scores{$locus}{"RAPTORX"}{"tmscore"}[$i]){
					my ($tmscore,$rmsd,$length,$seqid) = @{$alignment_scores{$locus}{"RAPTORX"}{"tmscore"}[$i]};
					$line .= "\t".$tmscore."\t".$rmsd."\t".$length;
				}
				else{
					$line .= "\t"x6;
				}
				if($alignment_scores{$locus}{"RAPTORX"}{"gdt"}[$i]){
					my ($gdt,$rmsd,$nalign) = @{$alignment_scores{$locus}{"RAPTORX"}{"gdt"}[$i]};
					$line .= "\t".$gdt."\t".$rmsd."\t".$nalign;
				}
				else{
					$line .= "\t"x6;
				}
				$line .= "\t".$description;
			}
			else{
				$line .= "\t"x7;
				$skip++;
			}
			if($skip >= 2){
				print("\n");
				next ANNOT;
			}
			print("$line\n");
		}
	}
}