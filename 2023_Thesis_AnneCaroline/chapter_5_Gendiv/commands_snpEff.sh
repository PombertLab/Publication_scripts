#!/usr/bin/env bash

GENDIV=/home/jpombert/Analyses/Encephalitozoon/GenDiv

###### FastANI #####
ls $GENDIV/FASTA/*.fasta > $GENDIV/FASTA/fasta.list
fastANI \
  --rl $GENDIV/FASTA/fasta.list \
  --ql $GENDIV/FASTA/fasta.list \
  -o $GENDIV/fastANI.txt

###### Mash #####
run_Mash.pl \
  -f $GENDIV/FASTA/*.fasta \
  -o $GENDIV/MASH

###### Mummer ######
# hellem
for y in {50451,50504,50604,Swiss}; do
	REF="$GENDIV/FASTA/$y.fasta";
	for x in {50451,50504,50604,Swiss}; do
		QUERY="$GENDIV/FASTA/$x.fasta";
		PREFIX="$x.vs.$y";
		nucmer \
		  --delta $GENDIV/Mummer/$PREFIX.delta \
		  $REF \
		  $QUERY;
		dnadiff \
		-d $GENDIV/Mummer/$PREFIX.delta \
		-p $GENDIV/Mummer/$PREFIX;
	done;
done

# intestinalis
for y in {50506,50507,50603,50651}; do
	REF="$GENDIV/FASTA/$y.fasta";
	for x in {50506,50507,50603,50651}; do
		QUERY=".$GENDIV/FASTA/$x.fasta";
		PREFIX="$x.vs.$y";
		nucmer \
		  --delta $GENDIV/Mummer/$PREFIX.delta \
		  $REF \
		  $QUERY;
		dnadiff \
		-d $GENDIV/Mummer/$PREFIX.delta \
		-p $GENDIV/Mummer/$PREFIX;
	done;
done

# cuniculi
for y in {50602,ECI,ECII,ECIII,ECII_CZ,ECIII_L,GB_M1}; do
    REF="$GENDIV/FASTA/$y.fasta";
	for x in {50602,ECI,ECII,ECIII,ECII_CZ,ECIII_L,GB_M1}; do
		QUERY="$GENDIV/FASTA/$x.fasta";
		PREFIX="$x.vs.$y";
		nucmer \
		   --delta $GENDIV/Mummer/$PREFIX.delta \
		   $REF \
		   $QUERY;
		dnadiff \
		-d $GENDIV/Mummer/$PREFIX.delta \
		-p $GENDIV/Mummer/$PREFIX;
	done;
done

# Ordospora
for y in {OC4,GBEP,NOV7,FISK}; do
    REF="$GENDIV/FASTA/$y.fasta";
	for x in {OC4,GBEP,NOV7,FISK}; do
		QUERY="$GENDIV/FASTA/$x.fasta";
		PREFIX="$x.vs.$y";
		nucmer \
		   --delta $GENDIV/Mummer/$PREFIX.delta \
		   $REF \
		   $QUERY;
		dnadiff \
		-d $GENDIV/Mummer/$PREFIX.delta \
		-p $GENDIV/Mummer/$PREFIX;
	done;
done

########################################### READ MAPPING #######################################################

######################## SSRG.pl 

SSRG.pl \
  -f $GENDIV/FASTA/*.fasta \
  -o $GENDIV/FASTQ/FASTQ_150 \
  -gzip \
  -r 150
 
SSRG.pl \
  -f $GENDIV/FASTA/*.fasta \
  -o $GENDIV/FASTQ/FASTQ_250 \
  -gzip \
  -r 250

######################## get_SNPs.pl 

######## SSRG 150 bp
# hellem
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50451,50504,50604,Swiss}.fasta \
  -pe1 $GENDIV/FASTQ/FASTQ_150/{50451,50504,50604,Swiss}.150.R1.fastq.gz \
  -pe2 $GENDIV/FASTQ/FASTQ_150/{50451,50504,50604,Swiss}.150.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/SSRG_150

# intestinalis
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50506,50507,50603,50651}.fasta \
  -pe1 $GENDIV/FASTQ/FASTQ_150/{50506,50507,50603,50651}.150.R1.fastq.gz \
  -pe2 $GENDIV/FASTQ/FASTQ_150/{50506,50507,50603,50651}.150.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/SSRG_150

# cuniculi
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50602,ECII_CZ,ECIII_L,GB_M1,ECI,ECII,ECIII}.fasta \
  -pe1 $GENDIV/FASTQ/FASTQ_150/{50602,ECII_CZ,ECIII_L,GB_M1,ECI,ECII,ECIII}.150.R1.fastq.gz \
  -pe2 $GENDIV/FASTQ/FASTQ_150/{50602,ECII_CZ,ECIII_L,GB_M1,ECI,ECII,ECIII}.150.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/SSRG_150

# Ordospora
get_SNPs.pl \
  -fa $GENDIV/FASTA/{OC4,GBEP,FISK,NOV7}.fasta \
  -pe1 $GENDIV/FASTQ/FASTQ_150/{OC4,GBEP,FISK,NOV7}.150.R1.fastq.gz \
  -pe2 $GENDIV/FASTQ/FASTQ_150/{OC4,GBEP,FISK,NOV7}.150.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/SSRG_150


######## SSRG 250 bp
# hellem
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50451,50504,50604,Swiss}.fasta \
  -pe1 $GENDIV/FASTQ/FASTQ_250/{50451,50504,50604,Swiss}.250.R1.fastq.gz \
  -pe2 $GENDIV/FASTQ/FASTQ_250/{50451,50504,50604,Swiss}.250.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/SSRG_250

# cuniculi
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50602,ECII_CZ,ECIII_L,GB_M1,ECI,ECII,ECIII}.fasta \
  -pe1 $GENDIV/FASTQ/FASTQ_250/{50602,ECII_CZ,ECIII_L,GB_M1,ECI,ECII,ECIII}.250.R1.fastq.gz \
  -pe2 $GENDIV/FASTQ/FASTQ_250/{50602,ECII_CZ,ECIII_L,GB_M1,ECI,ECII,ECIII}.250.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/SSRG_250

# intestinalis
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50506,50507,50603,50651}.fasta \
  -pe1 $GENDIV/FASTQ/FASTQ_250/{50506,50507,50603,50651}.250.R1.fastq.gz \
  -pe2 $GENDIV/FASTQ/FASTQ_250/{50506,50507,50603,50651}.250.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/SSRG_250

# Ordospora
get_SNPs.pl \
  -fa $GENDIV/FASTA/{OC4,GBEP,FISK,NOV7}.fasta \
  -pe1 $GENDIV/FASTQ/FASTQ_250/{OC4,GBEP,FISK,NOV7}.250.R1.fastq.gz \
  -pe2 $GENDIV/FASTQ/FASTQ_250/{OC4,GBEP,FISK,NOV7}.250.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/SSRG_250

######## SRA / Fastp Illumina data 
#
# Test --min-var-freq 0.2 vs 0.7? Can we simplify? No, lowering to 0.2 generates
# too many false positives. must do heterozigosity independently

FASTP=$GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP

# hellem
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50451,50604,Swiss,50504}.fasta \
  -pe1 $FASTP/{50451,50604,Swiss}.*.R1.fastq.gz \
  -pe2 $FASTP/{50451,50604,Swiss}.*.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA

# cuniculi
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50602,GB_M1,ECII_CZ,ECIII_L,ECI,ECII,ECIII}.fasta \
  -pe1 $FASTP/{50602,ECIII_L,ECI,ECII,ECIII}.*.R1.fastq.gz \
  -pe2 $FASTP/{50602,ECIII_L,ECI,ECII,ECIII}.*.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA

get_SNPs.pl \
  -fa $GENDIV/FASTA/{50602,GB_M1,ECII_CZ,ECIII_L,ECI,ECII,ECIII}.fasta \
  -fq $FASTP/ECII_CZ.*.SE.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA

# intestinalis
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50506,50507,50603,50651}.fasta \
  -pe1 $FASTP/{50506,50507,50603,50651}.*.R1.fastq.gz \
  -pe2 $FASTP/{50506,50507,50603,50651}.*.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  --min-avg-qual 30 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA

# Ordospora
get_SNPs.pl \
  -fa $GENDIV/FASTA/{OC4,GBEP,NOV7,FISK}.fasta \
  -pe1 $FASTP/{GBEP,NOV7,FISK}.*.R1.fastq.gz \
  -pe2 $FASTP/{GBEP,NOV7,FISK}.*.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  --min-avg-qual 30 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA

######## copy + rename VCFs for SnpEff
$GENDIV/Scripts/rename_vcf.pl \
  -v $GENDIV/RM/SSRG_150/minimap2.varscan2.VCFs/*.vcf \
  -o $GENDIV/snpEFF/VCFs_150/

$GENDIV/Scripts/rename_vcf.pl \
  -v $GENDIV/RM/SSRG_250/minimap2.varscan2.VCFs/*.vcf \
  -o $GENDIV/snpEFF/VCFs_250/

$GENDIV/Scripts/rename_vcf.pl \
  -v $GENDIV/RM/ILLUMINA/minimap2.varscan2.VCFs/*.vcf \
  -o $GENDIV/snpEFF/VCFs_ILLUMINA/

###### SnpEff  5.1d    2022-04-19

SNPEFF=$GENDIV/snpEFF

mkdir $SNPEFF/data
mkdir $SNPEFF/data/hellem50604
cd $SNPEFF/data/hellem50604

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/399/255/GCA_024399255.1_ASM2439925v1/GCA_024399255.1_ASM2439925v1_genomic.fna.gz
mv GCA_024399255.1_ASM2439925v1_genomic.fna.gz hellem50604.fa.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/399/255/GCA_024399255.1_ASM2439925v1/GCA_024399255.1_ASM2439925v1_genomic.gbff.gz
mv GCA_024399255.1_ASM2439925v1_genomic.gbff.gz genes.gbk.gz
gunzip genes.gbk.gz

cat genes.gbk | grep -P '^VERSION' | sed 's/VERSION     //'

# CP075147.1, CP075148.1, CP075149.1, CP075150.1, CP075151.1, CP075152.1, CP075153.1, CP075154.1, CP075155.1, CP075156.1, CP075157.1

## Add to $SNPEFF/snpEff.config

# Encephalitozoon hellem ATCC 50604
50604.genome : hellem
50604.chromosomes : CP075147.1, CP075148.1, CP075149.1, CP075150.1, CP075151.1, CP075152.1, CP075153.1, CP075154.1, CP075155.1, CP075156.1, CP075157.1

# Encephalitozoon hellem ATCC 50451
50451.genome : hellem
50451.chromosomes : CP119062.1, CP119063.1, CP119064.1, CP119065.1, CP119066.1, CP119067.1, CP119068.1, CP119069.1, CP119070.1, CP119071.1, CP119072.1

# Encephalitozoon hellem Swiss
Swiss.genome : hellem
Swiss.chromosomes : JAGEXW010000001.1, JAGEXW010000002.1, JAGEXW010000003.1, JAGEXW010000004.1, JAGEXW010000005.1, JAGEXW010000006.1, JAGEXW010000007.1, JAGEXW010000008.1, JAGEXW010000009.1, JAGEXW010000010.1, JAGEXW010000011.1, JAGEXW010000012.1, JAGEXW010000013.1, JAGEXW010000014.1, JAGEXW010000015.1, JAGEXW010000016.1, JAGEXW010000017.1, JAGEXW010000018.1, JAGEXW010000019.1, JAGEXW010000020.1, JAGEXW010000021.1, JAGEXW010000022.1, JAGEXW010000023.1, JAGEXW010000024.1, JAGEXW010000025.1, JAGEXW010000026.1, JAGEXW010000027.1, JAGEXW010000028.1, JAGEXW010000029.1, JAGEXW010000030.1, JAGEXW010000031.1, JAGEXW010000032.1

# Encephalitozoon hellem ATCC 50504
50504.genome : hellem
50504.chromosomes : NC_018468.1, NC_018469.1, NC_018470.1, NC_018471.1, NC_018472.1, NC_018473.1, NC_018474.1, NC_018475.1, NW_004057897.1, NW_004057898.1 NC_018479.1, NC_018480.1

# Encephalitozoon intestinalis ATCC 50506
50506.genome : intestinalis
50506.chromosomes : CP075158.1, CP075159.1, CP075160.1, CP075161.1, CP075162.1, CP075163.1, CP075164.1, CP075165.1, CP075166.1, CP075167.1, CP075168.1

# Encephalitozoon intestinalis ATCC 50507
50507.genome : intestinalis
50507.chromosomes : contig_01, contig_02, contig_03, contig_04, contig_05, contig_06, contig_07, contig_08, contig_09, contig_10, contig_11, contig_12, contig_13, contig_14, contig_15, contig_16, contig_17, contig_18, contig_19, contig_20, contig_21, contig_22, contig_23, contig_24, contig_25, contig_26, contig_27, contig_28

# Encephalitozoon intestinalis ATCC 50603
50603.genome : intestinalis
50603.chromosomes : contig_01, contig_02, contig_03, contig_04, contig_05, contig_06, contig_07, contig_08, contig_09, contig_10, contig_11, contig_12, contig_13, contig_14, contig_15, contig_16, contig_17, contig_18, contig_19, contig_20, contig_21, contig_22, contig_23, contig_24, contig_25, contig_26, contig_27, contig_28, contig_29, contig_30, contig_31, contig_32

# Encephalitozoon intestinalis ATCC 50651
50651.genome : intestinalis
50651.chromosomes : contig_01, contig_02, contig_03, contig_04, contig_05, contig_06, contig_07, contig_08, contig_09, contig_10, contig_11, contig_12, contig_13, contig_14, contig_15, contig_16, contig_17, contig_18, contig_19, contig_20, contig_21, contig_22, contig_23, contig_24, contig_25, contig_26, contig_27, contig_28

# Encephalitozoon cuniculi GB-M1
GB_M1.genome : cuniculi
GB_M1.chromosomes : NC_003242.2, NC_003229.1, NC_003230.1, NC_003231.1, NC_003232.1, NC_003233.1, NC_003234.1, NC_003235.1, NC_003238.2, NC_003236.1, NC_003237.2

# Encephalitozoon cuniculi ATCC 50602
50602.genome : cuniculi
50602.chromosomes : CP091441.1, CP091440.1, CP091439.1, CP091436.1, CP091438.1, CP091437.1, CP091435.1, CP091434.1, CP091433.1, CP091432.1, CP091431.1

# Encephalitozoon cuniculi ECIII_L
ECIII_L.genome : cuniculi
ECIII_L.chromosomes : LFTZ01000001.1, LFTZ01000002.1, LFTZ01000003.1, LFTZ01000004.1, LFTZ01000005.1, LFTZ01000006.1, LFTZ01000007.1, LFTZ01000008.1, LFTZ01000009.1, LFTZ01000010.1, LFTZ01000011.1, LFTZ01000012.1, LFTZ01000013.1, LFTZ01000014.1, LFTZ01000015.1

# Encephalitozoon cuniculi ECII_CZ
ECII_CZ.genome : cuniculi
ECII_CZ.chromosomes : KC513604.1, KC513605.1, KC513606.1, KC513607.1, KC513608.1, KC513609.1, KC513610.1, KC513611.1, KC513612.1, KC513613.1, KC513614.1, KC513615.1, KC513616.1, KC513617.1, KC513618.1, KC513619.1, KC513620.1, KC513621.1, KC513622.1, KC513623.1, KC513624.1, KC513625.1, KC513626.1, KC513627.1, KC513628.1, KC513629.1, KC513630.1, KC513631.1, KC513632.1, KC513633.1, KC513634.1, KC513635.1, KC513636.1, KC513637.1, KC513638.1, KC513639.1, KC513640.1, KC513641.1, KC513642.1, KC513643.1, KC513644.1, KC513645.1, KC513646.1, KC513647.1, KC513648.1, KC513649.1, KC513650.1, KC513651.1, KC513652.1, KC513653.1, KC513654.1, KC513655.1, KC513656.1, KC513657.1

# Encephalitozoon cuniculi ECI
ECI.genome : cuniculi
ECI.chromosomes : ECI_CH01, ECI_CH02, ECI_CH03, ECI_CH04, ECI_CH05, ECI_CH06, ECI_CH07, ECI_CH08, ECI_CH09a, ECI_CH09b, ECI_CH10, ECI_CH11

# Encephalitozoon cuniculi ECII
ECII.genome : cuniculi
ECII.chromosomes : ECII_CH01, ECII_CH02, ECII_CH03, ECII_CH04, ECII_CH05, ECII_CH06, ECII_CH07, ECII_CH08, ECII_CH09a, ECII_CH09b, ECII_CH09c,  ECII_CH09d, ECII_CH09e, ECII_CH10, ECII_CH11

# Encephalitozoon cuniculi ECIII
ECIII.genome : cuniculi
ECIII.chromosomes : ECIII_CH01, ECIII_CH02, ECIII_CH03, ECIII_CH04, ECIII_CH05, ECIII_CH06, ECIII_CH07, ECIII_CH08, ECIII_CH09a, ECIII_CH09b, ECIII_CH09c, ECIII_CH09d, ECIII_CH09e, ECIII_CH10, ECIII_CH11

# Ordospora colligata OC4
OC4.genome : ordospora
OC4.chromosomes : NW_014575407.1, NW_014575406.1, NW_014575405.1, NW_014575404.1, NW_014575403.1, NW_014575402.1, NW_014575401.1, NW_014575400.1, NW_014575399.1, NW_014575398.1, NW_014575397.1, NW_014575396.1, NW_014575395.1, NW_014575394.1, NW_014575393.1

# Ordospora colligata NOV7
NOV7.genome : ordospora
NOV7.chromosomes : PITF01000001.1, PITF01000002.1, PITF01000003.1, PITF01000004.1, PITF01000005.1, PITF01000006.1, PITF01000007.1, PITF01000008.1, PITF01000009.1, PITF01000010.1, PITF01000011.1, PITF01000012.1, PITF01000013.1, PITF01000014.1, PITF01000015.1, PITF01000016.1, PITF01000017.1, PITF01000018.1, PITF01000019.1, PITF01000020.1, PITF01000021.1

# Ordospora colligata FISK
FISK.genome : ordospora
FISK.chromosomes : PITH01000001.1, PITH01000002.1, PITH01000003.1, PITH01000004.1, PITH01000005.1, PITH01000006.1, PITH01000007.1, PITH01000008.1, PITH01000009.1, PITH01000010.1, PITH01000011.1, PITH01000012.1, PITH01000013.1, PITH01000014.1, PITH01000015.1, PITH01000016.1, PITH01000017.1, PITH01000018.1, PITH01000019.1, PITH01000020.1, PITH01000021.1, PITH01000022.1, PITH01000023.1, PITH01000024.1, PITH01000025.1, PITH01000026.1

# Ordospora colligata GBEP
GBEP.genome : ordospora
GBEP.chromosomes : PITG01000001.1, PITG01000002.1, PITG01000003.1, PITG01000004.1, PITG01000005.1, PITG01000006.1, PITG01000007.1, PITG01000008.1, PITG01000009.1, PITG01000010.1, PITG01000011.1, PITG01000012.1, PITG01000013.1, PITG01000014.1, PITG01000015.1, PITG01000016.1, PITG01000017.1, PITG01000018.1


## Microsporidiadb from gff HOWTO
wget https://microsporidiadb.org/common/downloads/release-62/EcuniculiEC2/fasta/data/MicrosporidiaDB-62_EcuniculiEC2_AnnotatedCDSs.fasta
mv MicrosporidiaDB-62_EcuniculiEC2_AnnotatedCDSs.fasta cds.fa
wget https://microsporidiadb.org/common/downloads/release-62/EcuniculiEC2/fasta/data/MicrosporidiaDB-62_EcuniculiEC2_AnnotatedProteins.fasta
mv MicrosporidiaDB-62_EcuniculiEC2_AnnotatedProteins.fasta protein.fa
sed 's/-p1//g' protein.fa > 1.fa
mv 1.fa  protein.fa
wget https://microsporidiadb.org/common/downloads/release-62/EcuniculiEC2/fasta/data/MicrosporidiaDB-62_EcuniculiEC2_Genome.fasta
mv MicrosporidiaDB-62_EcuniculiEC2_Genome.fasta sequences.fa
wget https://microsporidiadb.org/common/downloads/release-62/EcuniculiEC2/gff/data/MicrosporidiaDB-62_EcuniculiEC2.gff
mv MicrosporidiaDB-62_EcuniculiEC2.gff genes.gff

## build databases
java -jar $SNPEFF/snpEff.jar build -genbank 50604
java -jar $SNPEFF/snpEff.jar build -genbank 50451
java -jar $SNPEFF/snpEff.jar build -genbank Swiss
java -jar $SNPEFF/snpEff.jar build -genbank 50504
java -jar $SNPEFF/snpEff.jar build -genbank 50506
java -jar $SNPEFF/snpEff.jar build -genbank 50507
java -jar $SNPEFF/snpEff.jar build -genbank 50603
java -jar $SNPEFF/snpEff.jar build -genbank 50651
java -jar $SNPEFF/snpEff.jar build -genbank GB_M1
java -jar $SNPEFF/snpEff.jar build -genbank 50602
java -jar $SNPEFF/snpEff.jar build -genbank ECIII_L
java -jar $SNPEFF/snpEff.jar build -genbank ECII_CZ
java -jar $SNPEFF/snpEff.jar build -gff3 -v ECI
java -jar $SNPEFF/snpEff.jar build -gff3 -v ECII
java -jar $SNPEFF/snpEff.jar build -gff3 -v ECIII
java -jar $SNPEFF/snpEff.jar build -genbank OC4
java -jar $SNPEFF/snpEff.jar build -genbank NOV7
java -jar $SNPEFF/snpEff.jar build -genbank FISK
java -jar $SNPEFF/snpEff.jar build -genbank GBEP

## VCF name structure - hellem50451.vs.hellem50451.vcf
## Run snpEff with html
for y in {50451,50504,50604,Swiss}; do
	REF=$y;
	for x in {50451,50504,50604,Swiss}; do
		VCF="$x.vs.$REF.vcf";
		SUM="$x.vs.$REF.summary.html";
		OUT="$x.vs.$REF.annotations.vcf";
		java -Xmx8g -jar $SNPEFF/snpEff.jar \
		  -c $SNPEFF/snpEff.config \
		  -v \
		  $REF \
		  $VCF \
		  -s $SUM \
		  > $OUT
	done;
done

#####  Run snpEff with csv stats

SNPEFF=$GENDIV/snpEFF
for dir in {150,250,ILLUMINA}; do

	cd $SNPEFF/VCFs_$dir;

	# hellem
	for y in {50451,50504,50604,Swiss}; do
		REF="$y";
		for x in {50451,50504,50604,Swiss}; do
			VCF="$x.vs.$y.vcf";
			SUM="$x.vs.$y.summary.csv";
			OUT="$x.vs.$y.annotations.vcf";
			java -Xmx8g -jar $SNPEFF/snpEff.jar \
			  -c $SNPEFF/snpEff.config \
			  -v \
			  $REF \
			  $VCF \
			  -csvStats $SUM  \
			  > $OUT
		done;
	done

	# intestinalis
	for y in {50506,50507,50603,50651}; do
		REF="$y";
		for x in {50506,50507,50603,50651}; do
			VCF="$x.vs.$y.vcf";
			SUM="$x.vs.$y.summary.csv";
			OUT="$x.vs.$y.annotations.vcf";
			java -Xmx8g -jar $SNPEFF/snpEff.jar \
			  -c $SNPEFF/snpEff.config \
			  -v \
			  $REF \
			  $VCF \
			  -csvStats $SUM  \
			  > $OUT
		done;
	done

	# cuniculi
	for y in {50602,ECII_CZ,ECIII_L,GB_M1,ECI,ECII,ECIII}; do
		REF="$y";
		for x in {50602,ECII_CZ,ECIII_L,GB_M1,ECI,ECII,ECIII}; do
			VCF="$x.vs.$y.vcf";
			SUM="$x.vs.$y.summary.csv";
			OUT="$x.vs.$y.annotations.vcf";
			java -Xmx8g -jar $SNPEFF/snpEff.jar \
			  -c $SNPEFF/snpEff.config \
			  -v \
			  $REF \
			  $VCF \
			  -csvStats $SUM  \
			  > $OUT
		done;
	done

	# Ordospora
	for y in {GBEP,NOV7,FISK,OC4}; do
		REF="$y";
		for x in {GBEP,NOV7,FISK,OC4}; do
			VCF="$x.vs.$y.vcf";
			SUM="$x.vs.$y.summary.csv";
			OUT="$x.vs.$y.annotations.vcf";
			java -Xmx8g -jar $SNPEFF/snpEff.jar \
			  -c $SNPEFF/snpEff.config \
			  -v \
			  $REF \
			  $VCF \
			  -csvStats $SUM  \
			  > $OUT
		done;
	done

done

## make the table
$GENDIV/Scripts/make_table_SNP.pl \
  -list $GENDIV/Scripts/Table_list_SNP.txt \
  -acc $GENDIV/Scripts/accessions.tsv \
  -ani $GENDIV/fastANI.txt \
  -mash $GENDIV/MASH/Mash.mash \
  -mummer $GENDIV/Mummer/*.report \
  -s150 $GENDIV/RM/SSRG_150/minimap2.varscan2.stats/*.stats \
  -s250 $GENDIV/RM/SSRG_250/minimap2.varscan2.stats/*.stats \
  -sill $GENDIV/RM/ILLUMINA/minimap2.varscan2.stats/*.stats \
  -e150 $GENDIV/snpEFF/VCFs_150/*.csv \
  -e250 $GENDIV/snpEFF/VCFs_250/*.csv \
  -eill $GENDIV/snpEFF/VCFs_ILLUMINA/*.csv \
  -out $GENDIV/table_SNP.tsv

### Prototypes ...
java -Xmx8g -jar $SNPEFF/snpEff.jar \
  -c $SNPEFF/snpEff.config \
  -v \
  50506 \
  50506.vs.50506.vcf \
  -csvStats 50506.vs.50506.summary.csv  \
  > 50506.vs.50506.annotations.vcf
  