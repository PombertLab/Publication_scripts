#!/usr/bin/perl
## Pombert Lab, Illinois Tech, ACMdS, Feb 2022

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;

my $name = 'chrom_table.pl';
my $version = '0.3';
my $update = '2022-11-30';

my $usage = <<"USAGE";

NAME        $name
VERSION     $version
UPDATE      $update
SYNOPSIS    Build table of telomere, subtelomere and core nucleotide positions. \
            Calculate deltas of chromosome parts. \
            Counts rRNAs, tRNAs, CDS, core genes, subtelomere genes and total genes.

EXAMPLE     ./$name -fa eint_genome.fasta -gc eint_gc_shift.tsv -gb eint.gb -t eint_TideHunter.tsv -o out.tsv

OPTIONS
-f          .fasta file with all chromosomes
-t          .tsv file from TideHunder
-gc         .tsv file with GC shift position
-gb         Genbank (.gb) file
-o          output file name [Default: out.tsv]

USAGE

die "$usage\n" unless @ARGV;

my $fasta;
my $tide_hunter;
my $gc;
my $gb;
my $out = 'out.tsv';
my $db;

GetOptions(
    'f=s' => \$fasta,
    't=s' => \$tide_hunter,
    'gc=s' => \$gc,
    'gb=s' => \$gb,
    'o=s' => \$out,
);

###########################################################################################
###### DEFINING DATABASES
###########################################################################################

# Database of E. intestinalis, E. hellem and E. cuniculi chromosome locus_tags and accession numbers
my %L2C = (
    'GPK93_01' => 'CP075158',
    'GPK93_02' => 'CP075159',
    'GPK93_03' => 'CP075160',
    'GPK93_04' => 'CP075161',
    'GPK93_05' => 'CP075162',
    'GPK93_06' => 'CP075163',
    'GPK93_07' => 'CP075164',
    'GPK93_08' => 'CP075165',
    'GPK93_09' => 'CP075166',
    'GPK93_10' => 'CP075167',
    'GPK93_11' => 'CP075168',
    'GPU96_01' => 'CP075147',
    'GPU96_02' => 'CP075148',
    'GPU96_03' => 'CP075149',
    'GPU96_04' => 'CP075150',
    'GPU96_05' => 'CP075151',
    'GPU96_06' => 'CP075152',
    'GPU96_07' => 'CP075153',
    'GPU96_08' => 'CP075154',
    'GPU96_09' => 'CP075155',
    'GPU96_10' => 'CP075156',
    'GPU96_11' => 'CP075157',
    'JA071_11' => 'CP091441',
    'JA071_10' => 'CP091440',
    'JA071_09' => 'CP091439',
    'JA071_06' => 'CP091436',
    'JA071_08' => 'CP091438',
    'JA071_07' => 'CP091437',
    'JA071_05' => 'CP091435',
    'JA071_04' => 'CP091434',
    'JA071_03' => 'CP091433',
    'JA071_02' => 'CP091432',
    'JA071_01' => 'CP091431'
);

# Database of E. intestinalis, E. hellem and E. cuniculi chromosome accession numbers and numbers
my %C2C = (
    'CP075147' => 1,
    'CP075148' => 2,
    'CP075149' => 3,
    'CP075150' => 4,
    'CP075151' => 5,
    'CP075152' => 6,
    'CP075153' => 7,
    'CP075154' => 8,
    'CP075155' => 9,
    'CP075156' => 10,
    'CP075157' => 11,
    'CP091441' => 1,
    'CP091440' => 2,
    'CP091439' => 3,
    'CP091436' => 4,
    'CP091438' => 5,
    'CP091437' => 6,
    'CP091435' => 7,
    'CP091434' => 8,
    'CP091433' => 9,
    'CP091432' => 10,
    'CP091431' => 11,
    'CP075158' => 1,
    'CP075159' => 2,
    'CP075160' => 3,
    'CP075161' => 4,
    'CP075162' => 5,
    'CP075163' => 6,
    'CP075164' => 7,
    'CP075165' => 8,
    'CP075166' => 9,
    'CP075167' => 10,
    'CP075168' => 11
);

###########################################################################################
###### FASTA FILE
###########################################################################################

my %sizes;
my $chrom_num;
my $gen_size;

open FA, "<", "$fasta";

while(my $line = <FA>){
    chomp $line;

    # Grabbing chromosme number
    if ($line =~ /^>(\w+)/){
        $chrom_num = $1;
    }
    # Assigning DNA sequence to each chromosome
    else {
        $sizes{$chrom_num}{'DNA'} .= $line;
    }
}

# Calculating the length of each chromosome and total genome from FASTA
for (sort (keys %sizes)){
    my $chrom = $_;
    my $DNA = $sizes{$chrom}{'DNA'};
    my $chrom_size = length($DNA);
    $sizes{$chrom}{'Length'} = $chrom_size;
    $gen_size += $sizes{$chrom}{'Length'};
}

###########################################################################################
###### GC SHIFT FILE
###########################################################################################

my %core;
my %subtelomere;

open GC, "<", "$gc";

while(my $line = <GC>){
    chomp $line;

    if ($line =~ /^#/){
        next;
    }
    else{
        # Spliting columns in tsv per tab
        my @columns = split("\t", $line);

        my $chrom = $columns[0];
        my $core_start = $columns[1];
        my $core_end = $columns[2];
        
        # Defining core start/end positions
        $core{$chrom}{'Core_start'} = $core_start;
        $core{$chrom}{'Core_end'} = $core_end;

        # Calculating subtelomeres start and end positions
        # Send is the end of first subtelomere | Estart it the start of the last subtelomere
        my $Send_subtel = $core_start -1;
        my $Estart_subtel = $core_end +1;

        $subtelomere{$chrom}{'Send'} = $Send_subtel;
        $subtelomere{$chrom}{'Estart'} = $Estart_subtel;

    }
}

###########################################################################################
###### GENBANK FILE
###########################################################################################

# Grabbing feature count and rRNA subtunits locations from GB file
my %rRNA;
my %features;
my $chrom;
my $FEAT;
my $type;
my $note_type;

my $gene_start;
my $gene_end;

# Variables for counting features for the genome
my $rRNA_count = 0;
my $tRNA_count = 0;
my $CDS_count = 0;
my $core_count = 0;
my $subtel_count = 0;
my $gene_count = 0;

# Variables for counting variables per chromosome
my $rRNA_chrm;
my $tRNA_chrm;
my $CDS_chrm;
my $gene_core;
my $gene_subtel;
my $gene_chrm;

open GB, "<", "$gb";

while (my $line = <GB>){
    chomp $line;
    
    # Grabing chromosome number and restarting variables per chromosome
    if($line =~ /^LOCUS\s+(\w+)/){
        $chrom = $1;
        $gene_chrm = 0;
        $gene_subtel = 0;
        $gene_core = 0;
        $rRNA_chrm = 0;
        $tRNA_chrm = 0;
        $CDS_chrm = 0;
    }
    # Working of finding features 
    elsif ($line =~ /^FEATURES/){
        $FEAT = 1;
        next;
    }
    #clean FEAT at the end of each chromosome
    elsif ($line =~ /ORIGIN/){
        $FEAT = undef;
    }
    if ($FEAT){
        my @rRNA_start;
        my @rRNA_end;

        # Determining which feature per gene
        if ($line =~ /^     (CDS|gene|source|rRNA|mRNA|tRNA)\s+.*?(\d+)\.\.>*(\d+)/){
            $type = $1;
            $gene_start = $2;
            $gene_end = $3;

            # Counting how many genes (all features included) per chromosome and genome
            if ($type eq 'gene'){
                $gene_count++;
                $gene_chrm++;
                
                # Counting genes in core and subtelomeres
                # First, genes which start and end positions are within the core boundaries
                if ($gene_start >= $core{$chrom}{'Core_start'} && $gene_end <= $core{$chrom}{'Core_end'}){
                    $core_count++;
                    $gene_core++;
                }
                # Second, genes with start position before the end of subtelomere 1 and after start of subtelomere 2
                elsif ($gene_start <= $subtelomere{$chrom}{'Estart'} or $gene_start >= $subtelomere{$chrom}{'Send'}){
                    $subtel_count++;
                    $gene_subtel++;
                }
            }
            # Counting how many tRNAs per chromosome and genome
            if($type eq 'tRNA'){
                $tRNA_count++;
                $tRNA_chrm++;
            }
            # Counting mRNAs or CDS per chromosome and genome
            if ($type eq 'CDS'){
                $CDS_count++;
                $CDS_chrm++;
            }
            # Counting rRNA per chromosome and genome
            if ($type eq 'rRNA'){
                $rRNA_count++;
                $rRNA_chrm++;
            }
            # Populating features database for number of features in general
            $features{$chrom}{'rRNA'} = $rRNA_chrm;
            $features{$chrom}{'tRNA'} = $tRNA_chrm;
            $features{$chrom}{'CDS'} = $CDS_chrm;
            $features{$chrom}{'core'} = $gene_core;
            $features{$chrom}{'subtel'} = $gene_subtel;
            $features{$chrom}{'gene'} = $gene_chrm;
        }
        # Working on getting position of rRNA subunits
        if ($type eq 'rRNA'){
            if ($line =~ /^\s+\/(\w+)/){
                $note_type = $1;

                if ($note_type eq 'product' && $line =~ /^\s+\/product=\"(small|large)\ssubunit.*\"/){
                    # Adding start and end position of rRNA subiunits per chromosome
                    push ( @{ $rRNA { $chrom }{'rRNA_start'} }, "$gene_start");
                    push ( @{ $rRNA { $chrom }{'rRNA_end'} }, "$gene_end");
                }
            }
        }
    }
} 

# Defining the start and end of the rRNA subunit block
for my $chr (keys %rRNA){
    # Working on start of rRNA subunits block
    my @pos_start = @{ $rRNA{$chr}{'rRNA_start'} };
    my $rRNA_start_1 = $pos_start[0];
    my $rRNA_start_2 = $pos_start[2];

    # Working on end of rRNA subunits block
    my @pos_end = @{ $rRNA{$chr}{'rRNA_end'} };
    my $rRNA_end_1 = $pos_start[1];
    my $rRNA_end_2 = $pos_start[3];

    # Populating database for later printing
    $rRNA{$chr}{'rRNA_start_1'} = $rRNA_start_1;
    $rRNA{$chr}{'rRNA_start_2'} = $rRNA_start_2;
    $rRNA{$chr}{'rRNA_end_1'} = $rRNA_end_1;
    $rRNA{$chr}{'rRNA_end_2'} = $rRNA_end_2;

    ## Defining subtelomere edges
    $subtelomere{$chr}{'Sstart'} = $rRNA_end_1 + 1;
    $subtelomere{$chr}{'Eend'} = $rRNA_start_2 - 1 ;

}

# Using the selected database to translate the locus tag to chromosome accession number
for my $keys (sort (keys %L2C)){

    my $chrom_accession = $L2C{$keys};
    my $chromie = $C2C{$chrom_accession};

    # Chaging locus tag for chromosome accession number in core db
    if (exists $core{$keys}){
        $core{$chrom_accession} = $core{$keys};
    }

    # Changing locus tag for chromosome accession number in subtelomere db
    if (exists $subtelomere{$keys}){
        $subtelomere{$chrom_accession} = $subtelomere{$keys};
    }
}

###########################################################################################
###### TIDE HUNTER FILE
###########################################################################################

my %chromosomes;

open TH, "<", "$tide_hunter";

while (my $line = <TH>){
    chomp $line;

    # Grabing important info from the TideHunter output file 
    my @columns = split("\t", $line);

    my $chromosome = $columns[0];
    my ($digit) = $columns[1] =~ /rep(\d+)/;
    $digit = sprintf ("%02d", $digit);
    my $frep = 'rep'."$digit"; 
    my $start = $columns[4];
    my $end = $columns[5];
    
    # Setting database: chromosome name, repat number attached to start and end of repeat in the chromosome
    $chromosomes{$chromosome}{$frep}{'start'} = $start;
    $chromosomes{$chromosome}{$frep}{'end'} = $end;

}

# Defining the start and end of subtelomere1 and subtelomere2, respectively
for my $chrom (sort(keys %chromosomes)){

    my @reps = sort (keys %{$chromosomes{$chrom}});
    
    # Grabing first and last repeat of the chromosome (telomeric repeats)
    my $first_rep = $reps[0];
    my $last_rep = $reps[-1];

    # Defining start and end of first telomere
    my $end_tel_1 = $chromosomes{$chrom}{$first_rep}{'end'};
    
    # Defining the start and end of last telomere
    my $start_tel_2 = $chromosomes{$chrom}{$last_rep}{'start'};
}

###########################################################################################
###### OUTPUT FILE
###########################################################################################

open TSV, ">", "$out";

# Defining headers
my @header = (
    'Chromosome',
    'Accession',
    'Length',
    'Tel_start',
    'Tel_end',
    'rRNA_start_1',
    'rRNA_end_1',
    'Subtel_start1',
    'Subtel_end1',
    'Core_start',
    'Core_end',
    'Subtel_start2',
    'Subtel_end2',
    'rRNA_start_2',
    'rRNA_end_2',
    'Tel_start2',
    'Tel_end2',
    'Delta_tel1',
    'Delta_subtel1',
    'Delta_core',
    'Delta_subtel2',
    'Delta_tel2',
    'rRNA_count',
    'tRNA_count',
    'CDS_count',
    'Core_genes',
    'Subtel_genes',
    'Total_genes'
);

# Printing headers to output
print TSV "#";
foreach my $num (0..($#header-1)){
    print TSV "$header[$num]\t";
}
print TSV "$header[-1]\n";

# Sorting database by value (chromosome #) rather than key (chromosome accession)
my @chr = sort { $C2C{$a} <=> $C2C{$b} } keys %C2C;

# Calculating the length of each chromosome part
my %deltas;
while (my $chr = shift @chr){

    # Formating chromosome numbers to be in correct order
    my $val = sprintf("%02d", $C2C{$chr});
    
    # Cheking if the database was populated
    if (exists ($chromosomes{$chr})){

        # MAking array of reps to use just first and last ones
        my @keys =  sort ( keys %{$chromosomes{$chr}} );

        # Calculating size (Deltas) of chromosome parts
        $deltas{$chr}{'Delta_tel1'} = $rRNA{$chr}{'rRNA_start_1'} - $chromosomes{$chr}{$keys[0]}{'start'};
        $deltas{$chr}{'Delta_tel2'} = $chromosomes{$chr}{$keys[-1]}{'end'} - $rRNA{$chr}{'rRNA_end_2'};
        $deltas{$chr}{'Delta_subtel1'} = $subtelomere{$chr}{'Send'} - $rRNA{$chr}{'rRNA_start_1'};
        $deltas{$chr}{'Delta_subtel2'} = $rRNA{$chr}{'rRNA_end_2'} - $subtelomere{$chr}{'Estart'};
        $deltas{$chr}{'Delta_core'} = $core{$chr}{'Core_end'} - $core{$chr}{'Core_start'}; 
        
        # Printing data in out file
        # Chromosomes accession and length
        print TSV "$val\t";                                     # chromosome #
        print TSV "$chr\t";                                     # chromosome accession number
        print TSV "$sizes{$chr}{'Length'}\t";                   # Chromosome length from FASTA
        # Chromosomes partitioning
        print TSV "$chromosomes{$chr}{$keys[0]}{'start'}\t";    # start telomeric repeat 1
        print TSV "$chromosomes{$chr}{$keys[0]}{'end'}\t";      # end telomeric repeat 1
        print TSV "$rRNA{$chr}{'rRNA_start_1'}\t";              # start of rRNA subunits 1
        print TSV "$rRNA{$chr}{'rRNA_end_1'}\t";                # end of rRNA subunits 1
        print TSV "$subtelomere{$chr}{'Sstart'}\t";             # start subtelomere 1
        print TSV "$subtelomere{$chr}{'Send'}\t";               # end subtelomere 1
        print TSV "$core{$chr}{'Core_start'}\t";                # start of core
        print TSV "$core{$chr}{'Core_end'}\t";                  # end of core
        print TSV "$subtelomere{$chr}{'Estart'}\t";             # start subtelomere 2
        print TSV "$subtelomere{$chr}{'Eend'}\t";               # end subtelomere 2
        print TSV "$rRNA{$chr}{'rRNA_start_2'}\t";              # start of rRNA subunits 2
        print TSV "$rRNA{$chr}{'rRNA_end_2'}\t";                # end of rRNA subunits 2
        print TSV "$chromosomes{$chr}{$keys[-1]}{'start'}\t";   # start telomere 2
        print TSV "$chromosomes{$chr}{$keys[-1]}{'end'}\t";     # end telomere 2 
        # Chromosomes deltas (length of each part)
        print TSV "$deltas{$chr}{'Delta_tel1'}\t";              # delta telomere 1
        print TSV "$deltas{$chr}{'Delta_subtel1'}\t";           # delta subtelomere 1
        print TSV "$deltas{$chr}{'Delta_core'}\t";              # delta core
        print TSV "$deltas{$chr}{'Delta_subtel2'}\t";           # delta subtelomere 2
        print TSV "$deltas{$chr}{'Delta_tel2'}\t";              # delta telomere 2
        # Number of features per chromosome
        print TSV "$features{$chr}{'rRNA'}\t";                  # rRNA count per chromosome
        print TSV "$features{$chr}{'tRNA'}\t";                  # tRNA count per chromosome
        print TSV "$features{$chr}{'CDS'}\t";                   # CDS count per chromosome
        print TSV "$features{$chr}{'core'}\t";                  # Gene (ORF) count in the core per chromosome 
        print TSV "$features{$chr}{'subtel'}\t";                # Gene (ORF) count in the subtelomeres per chromosome
        print TSV "$features{$chr}{'gene'}\n";                  # Gene (ORF) count per chromosome
    }
}
# Printing length of genome and total count of features
print TSV "\tTotal\t$gen_size\t";                               # Full size of genome from summing chromosome lengths
print TSV "\t"x19;                                              # Printing TABS
print TSV "$rRNA_count\t$tRNA_count\t$CDS_count\t";             # rRNA, tRNA and CDS counts for the genome
print TSV "$core_count\t$subtel_count\t$gene_count\t";          # Core, subtelomere and gene counts for the genome