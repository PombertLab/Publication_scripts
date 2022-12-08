#!/usr/bin/perl
## Pombert Lab, Illinois Tech, ACMdS, 2022

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $name = 'get_sub_proteins.pl';
my $version = '0.2';
my $update = '2022-07-04';

my $usage = <<"USAGE";

NAME            $name
VERSION         $version
UPDATE          $update
SYNOPSIS        Get CDS with function located in the subteloemric region of chromosomes

EXAMPLE         ./$name -c chromosome.tsv -g example.gb -o out.tsv

OPTIONS
-c (--chrom)    Input chromosome components table generated with chrom_table.pl
-g (--gb)       Input Genbank .gb file
-o (--output)   Output file name [Default: out.tsv]
-no_hp          No hypotetical proteins in the output file
-no_duf         No DUF domain proteins in output file

USAGE

die "$usage\n" unless @ARGV;

my $chrom;
my $gb;
my $out = 'out.tsv';
my $no_hp;
my $no_duf;

GetOptions(
    'c|chrom=s' => \$chrom,
    'g|gb=s' => \$gb,
    'o|output=s' => \$out,
    'no_hp' => \$no_hp,
    'no_duf' => \$no_duf
);

open CHR, "<", "$chrom";

my %pos;
while (my $line = <CHR>){
    chomp $line;

    if ($line =~ /^#/){ next; }
    else{
        my @columns = split("\t", $line);

        # Grabing important info about chromosome structure
        my $chrm = $columns[0];
        my $chrm_acc = $columns[1];
        my $start_sub1 = $columns[6];
        my $end_sub1 = $columns[7];
        my $start_sub2 = $columns[10];
        my $end_sub2 = $columns[11];

        # Populating the database of subtelomere start and end posi for each chromosome
        $pos{$chrm_acc}{$chrm}{'start_1'} = $start_sub1;
        $pos{$chrm_acc}{$chrm}{'start_2'} = $start_sub2;
        $pos{$chrm_acc}{$chrm}{'end_1'} = $end_sub1;
        $pos{$chrm_acc}{$chrm}{'end_2'} = $end_sub2;
    }
}

open GB, "<", "$gb";

my %cds;
my $FEAT;
my $chrom_locus;
my $lc;
my $cds_start;
my $cds_end;
my $product;
my $type;
my $note_type;

while (my $line = <GB>){
    chomp $line;

    # Grabing Locus of cromosomes
    if ($line =~ /^LOCUS\s+(\w+)\s+/){
        $chrom_locus = $1;
        # $cds{$chom_locus};
    }

    ## Making flag to look for info only between feature and origin
    if ($line =~ /^FEATURES/){
        $FEAT = 1;
        next;
    }
    elsif ($line =~ /ORIGIN/){
        $FEAT = undef;
    }
    if ($FEAT){
        
        # Grabing the start and end position of the CDS
        if ($line =~ /^     (CDS|gene|source|rRNA|mRNA|tRNA)\s+.*?(\d+)\.\.(\d+)/){
            $type = $1;
            $cds_start = $2;
            $cds_end = $3;
        }

        # Grabing info only from CDS
        if ($type eq 'CDS'){
            
            # Identifying type of note, if it is locus_tag, product_translation, etc
            if ($line =~ /^\s+\/(\w+)/){
                $note_type = $1;

                # Grabing locus tag from the CDS
                if ($note_type eq 'locus_tag' && $line =~ /^\s+\/locus_tag=\"(\w+)\"/){
                    $lc = $1;
                    # Assigning start and end position of CDS
                    $cds{$lc}{$chrom_locus}{'start'} = $cds_start;
                    $cds{$lc}{$chrom_locus}{'end'} = $cds_end;
                }
                
                # Working on product name, it can be multi-line
                elsif ($note_type eq 'product' && $line =~ /^\s+\/product=\"(.*)/){
                    # $note_type = $1;
                    $product = $1;
                    $product =~ s/\"//;  # Removing quotes from product names that are in a single line
                    $cds{$lc}{$chrom_locus}{'product'} = $product;
                }
                # Working on product note type and lines that end with " (product name in multi-line)
                else{
                    if ($note_type eq 'product' && $line =~ /\"$/){
                        my $no_character = $line;
                        # Removing non-wanted characteres from line
                        $no_character =~ s/\s+/ /;
                        $no_character =~ s/"//;
                        # Appending second part of product name to the hash
                        $cds{$lc}{$chrom_locus}{'product'} .= $no_character;
                    }
                }
            }
        }
    }
}

## Working on products located only in subtelomere locations
open OUT, ">", "$out";
my %sub;

# Printing header on out file
print OUT "#Region\tChromosome\tLocus_tag\tStart_pos\tEnd_pos\tProduct\n";

# Sorting %pos by accession number
foreach my $acc (sort (keys %pos)){
    
    # Accessing the locus_tags in the hash
    my @chrom = (keys %{ $pos{$acc} });

    # Accessing position of subtelomeres in each chromosome
    while ( my $key = shift@chrom){
        my $start1 = $pos{$acc}{$key}{'start_1'};
        my $end1 = $pos{$acc}{$key}{'end_1'};

        my $start2 = $pos{$acc}{$key}{'start_2'};
        my $end2 =  $pos{$acc}{$key}{'end_2'};

        # Accessing protein start and end from GB file
        foreach my $loci (sort(keys %cds)){
            
            if (exists $cds{$loci}{$acc}){
                
                my $prot_start = $cds{$loci}{$acc}{'start'};
                my $prot_end = $cds{$loci}{$acc}{'end'};
                my $cds_product = $cds{$loci}{$acc}{'product'};

                ## Removing hypothetical and DUF proteins from the output file
                if ($no_hp && $no_duf){
                    # Grabbing CDS withing the subtelomere 1 region
                    if ($start1 <= $prot_start && $prot_start <= $end1) {
                        
                        # Getting rid of the hipothetical proteins in the output file
                        if ($cds_product eq 'hypothetical protein' || $cds_product =~ /^DUF/ ){next;}
                        
                        # Printing out only CDS with annotation
                        else{
                            print OUT "Sub1\t$acc\t$loci\t$prot_start\t$prot_end\t$cds_product\n";
                        }
                    }
                    # Grabing CDS within subtelomere 2 region
                    elsif ($start2 <= $prot_start && $prot_start <= $end2){
                        
                        # Getting rid of hypothetical proteins in the output file
                        if ($cds_product eq 'hypothetical protein' || $cds_product =~ /^DUF/){next;}
                        
                        # printing out only CDS with annotation
                        else{
                            print OUT "Sub2\t$acc\t$loci\t$prot_start\t$prot_end\t$cds_product\n";
                        }
                    }
                }
                # Removing only hypothetical proteins
                elsif ($no_hp){
                    # Grabbing CDS withing the subtelomere 1 region
                    if ($start1 <= $prot_start && $prot_start <= $end1) {
                        
                        # Getting rid of the hipothetical proteins in the output file
                        if ($cds_product eq 'hypothetical protein'){next;}
                        
                        # Printing out only CDS with annotation
                        else{
                            print OUT "Sub1\t$acc\t$loci\t$prot_start\t$prot_end\t$cds_product\n";
                        }
                    }
                    # Grabing CDS within subtelomere 2 region
                    elsif ($start2 <= $prot_start && $prot_start <= $end2){
                        
                        # Getting rid of hypothetical proteins in the output file
                        if ($cds_product eq 'hypothetical protein'){next;}
                        
                        # printing out only CDS with annotation
                        else{
                            print OUT "Sub2\t$acc\t$loci\t$prot_start\t$prot_end\t$cds_product\n";
                        }
                    }
                }
                # Removing only DUF proteins
                elsif ($no_duf){
                    if ($start1 <= $prot_start && $prot_start <= $end1) {
                            
                            # Getting rid of the hipothetical proteins in the output file
                            if ($cds_product =~ /^DUF/){next;}
                            
                            # Printing out only CDS with annotation
                            else{
                                print OUT "Sub1\t$acc\t$loci\t$prot_start\t$prot_end\t$cds_product\n";
                            }
                        }
                    # Grabing CDS within subtelomere 2 region
                    elsif ($start2 <= $prot_start && $prot_start <= $end2){
                        
                        # Getting rid of hypothetical proteins in the output file
                        if ($cds_product =~ /^DUF/){next;}
                        
                        # printing out only CDS with annotation
                        else{
                            print OUT "Sub2\t$acc\t$loci\t$prot_start\t$prot_end\t$cds_product\n";
                        }
                    }
                }
                # Printing out all CDS in the out file
                else {
                    # For CDS in subtelomere region 1
                    if ($start1 <= $prot_start && $prot_start <= $end1) {
                        print OUT "Sub1\t$acc\t$loci\t$prot_start\t$prot_end\t$cds_product\n";
                    }
                    # For CDS in subtelomere region 2
                    elsif ($start2 <= $prot_start && $prot_start <= $end2){
                        print OUT "Sub2\t$acc\t$loci\t$prot_start\t$prot_end\t$cds_product\n";
                    }
                }
            }
        }
    }
}

