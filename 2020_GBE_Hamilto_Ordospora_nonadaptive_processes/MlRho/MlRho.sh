### on MiniMe
/media/Data_3/Karen/MlRho/

### Data from NCBI
mkdir Hamiltosporidium_magnivora_BEOM2
cd Hamiltosporidium_magnivora_BEOM2/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/325/065/GCA_004325065.1_BEOM2_v1/GCA_004325065.1_BEOM2_v1_genomic.fna.gz
gunzip GCA_004325065.1_BEOM2_v1_genomic.fna.gz

mkdir Hamiltosporidium_magnivora_ILBN2
cd Hamiltosporidium_magnivora_ILBN2/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/325/035/GCA_004325035.1_ASM432503v1/GCA_004325035.1_ASM432503v1_genomic.fna.gz
gunzip GCA_004325035.1_ASM432503v1_genomic.fna.gz

mkdir Hamiltosporidium_tvaerminnensis_FIOER33
cd Hamiltosporidium_tvaerminnensis_FIOER33
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/325/045/GCA_004325045.1_FIOER33_v1/GCA_004325045.1_FIOER33_v1_genomic.fna.gz
gunzip GCA_004325045.1_FIOER33_v1_genomic.fna.gz

mkdir Hamiltosporidium_tvaerminnensis_ILG3
cd Hamiltosporidium_tvaerminnensis_ILG3/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/325/075/GCA_004325075.1_ILG3_v1/GCA_004325075.1_ILG3_v1_genomic.fna.gz
gunzip GCA_004325075.1_ILG3_v1_genomic.fna.gz

##### BEOM2 #####
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
BEOM2_R1.truseq.cutadapt.fastq BEOM2_R2.truseq.cutadapt.fastq \
BEOM2_R1.PE.trimmo.fastq BEOM2_R1.SE.trimmo.fastq \
BEOM2_R2.PE.trimmo.fastq BEOM2_R2.SE.trimmo.fastq \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:100

cat BEOM2_R1.PE.trimmo.fastq BEOM2_R1.SE.trimmo.fastq BEOM2_R2.PE.trimmo.fastq BEOM2_R2.SE.trimmo.fastq > BEOM2.trimmo.fastq

get_SNPs.pl -fa GCA_004325065.1_BEOM2_v1_genomic.fna \
-fq BEOM2.trimmo.fastq \
-mapper minimap2 -preset sr \
-rmo -bam
samtools index minimap2.BAM/BEOM2.trimmo.fastq.GCA_004325065.1_BEOM2_v1_genomic.fna.minimap2.bam ### Average SE depth 16.24

samtools view -b minimap2.BAM/BEOM2.trimmo.fastq.GCA_004325065.1_BEOM2_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I > theta.BEOM2.c4.SE.txt; mlRho -m 1000 -M 1005 >> theta.BEOM2.c4.SE.txt

samtools view -b minimap2.BAM/BEOM2.trimmo.fastq.GCA_004325065.1_BEOM2_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 8
mlRho -M 0 -I > theta.BEOM2.c8.SE.txt; mlRho -m 1000 -M 1005 >> theta.BEOM2.c8.SE.txt

samtools view -b minimap2.BAM/BEOM2.trimmo.fastq.GCA_004325065.1_BEOM2_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 16
mlRho -M 0 -I > theta.BEOM2.c16.SE.txt; mlRho -m 1000 -M 1005 >> theta.BEOM2.c16.SE.txt

##### ILBN2 #####
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
IL_BN_2_ATTCCT_S2_R1_001.fastq.gz IL_BN_2_ATTCCT_S2_R2_001.fastq.gz \
ILBN2_R1.PE.trimmo.fastq ILBN2_R1.SE.trimmo.fastq \
ILBN2_R2.PE.trimmo.fastq ILBN2_R2.SE.trimmo.fastq \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:100

cat ILBN2_R1.PE.trimmo.fastq ILBN2_R1.SE.trimmo.fastq ILBN2_R2.PE.trimmo.fastq ILBN2_R2.SE.trimmo.fastq > ILBN2.trimmo.fastq

get_SNPs.pl -fa GCA_004325035.1_ASM432503v1_genomic.fna \
-fq ILBN2.trimmo.fastq \
-mapper minimap2 -preset sr \
-bam -rmo
samtools index minimap2.BAM/ILBN2.trimmo.fastq.GCA_004325035.1_ASM432503v1_genomic.fna.minimap2.bam ### Average SE depth 11.76

samtools view -b minimap2.BAM/ILBN2.trimmo.fastq.GCA_004325035.1_ASM432503v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I  > theta.ILBN2.c4.SE.txt; mlRho -m 1000 -M 1005  >> theta.ILBN2.c4.SE.txt

samtools view -b minimap2.BAM/ILBN2.trimmo.fastq.GCA_004325035.1_ASM432503v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 8
mlRho -M 0 -I  > theta.ILBN2.c8.SE.txt; mlRho -m 1000 -M 1005  >> theta.ILBN2.c8.SE.txt

samtools view -b minimap2.BAM/ILBN2.trimmo.fastq.GCA_004325035.1_ASM432503v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 16
mlRho -M 0 -I  > theta.ILBN2.c16.SE.txt; mlRho -m 1000 -M 1005  >> theta.ILBN2.c16.SE.txt

##### FIOER33 #####
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
FIOER33_ACAGT_L001_R1_001_Hamiltosporidiu.fastq.gz FIOER33_ACAGT_L001_R1_001_Hamiltosporidiu.fastq.gz \
FIOER33_R1.PE.trimmo.fastq FIOER33_R1.SE.trimmo.fastq \
FIOER33_R2.PE.trimmo.fastq FIOER33_R2.SE.trimmo.fastq \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:100

cat FIOER33_R1.PE.trimmo.fastq FIOER33_R1.SE.trimmo.fastq FIOER33_R2.PE.trimmo.fastq FIOER33_R2.SE.trimmo.fastq > FIOER33.fastq

get_SNPs.pl -fa GCA_004325045.1_FIOER33_v1_genomic.fna \
-fq FIOER33.fastq \
-mapper minimap2 -preset sr \
-bam -rmo
samtools index minimap2.BAM/FIOER33.fastq.GCA_004325045.1_FIOER33_v1_genomic.fna.minimap2.bam ## Average SE depth 308.04X

samtools view -b minimap2.BAM/FIOER33.fastq.GCA_004325045.1_FIOER33_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I  > theta.FIOER33.c4.SE.txt; mlRho -m 1000 -M 1005  >> theta.FIOER33.c4.SE.txt

samtools view -b minimap2.BAM/FIOER33.fastq.GCA_004325045.1_FIOER33_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 8
mlRho -M 0 -I  > theta.FIOER33.c8.SE.txt; mlRho -m 1000 -M 1005  >> theta.FIOER33.c8.SE.txt

samtools view -b minimap2.BAM/FIOER33.fastq.GCA_004325045.1_FIOER33_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 16
mlRho -M 0 -I  > theta.FIOER33.c16.SE.txt; mlRho -m 1000 -M 1005  >> theta.FIOER33.c16.SE.txt

##### ILG3 #####
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
IL-G-3_R1.fastq IL-G-3_R2.fastq \
ILG3_R1.PE.trimmo.fastq ILG3_R1.SE.trimmo.fastq \
ILG3_R2.PE.trimmo.fastq ILG3_R2.SE.trimmo.fastq \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:100

# SE mode; doesn't work well, not sure why
cat ILG3_R1.PE.trimmo.fastq ILG3_R1.SE.trimmo.fastq ILG3_R2.PE.trimmo.fastq ILG3_R2.SE.trimmo.fastq > ILG3.trimmo.fastq

get_SNPs.pl -fa GCA_004325075.1_ILG3_v1_genomic.fna \
-fq ILG3.trimmo.fastq \
-mapper minimap2 -preset sr \
-bam -rmo
samtools index minimap2.BAM/ILG3.trimmo.fastq.GCA_004325075.1_ILG3_v1_genomic.fna.minimap2.bam ## Average SE depth 421.05X

samtools view -b minimap2.BAM/ILG3.trimmo.fastq.GCA_004325075.1_ILG3_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I  > theta.ILG3.c4.SE.txt; mlRho -m 1000 -M 1005  >> theta.ILG3.c4.SE.txt

samtools view -b minimap2.BAM/ILG3.trimmo.fastq.GCA_004325075.1_ILG3_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 8
mlRho -M 0 -I  > theta.ILG3.c8.SE.txt; mlRho -m 1000 -M 1005  >> theta.ILG3.c8.SE.txt

samtools view -b minimap2.BAM/ILG3.trimmo.fastq.GCA_004325075.1_ILG3_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 16
mlRho -M 0 -I  > theta.ILG3.c16.SE.txt; mlRho -m 1000 -M 1005  >> theta.ILG3.c16.SE.txt

# PE mode - 36X depth; works great at 36X, struggles at > 300X
#!/usr/bin/bash
split -l 24000000 \
--numeric-suffixes=01 \
--additional-suffix=.fastq \
IL-G-3_R1.fastq  \
R1_
split -l 24000000 \
--numeric-suffixes=01 \
--additional-suffix=.fastq \
IL-G-3_R2.fastq  \
R2_ 

get_SNPs.pl -fa GCA_004325075.1_ILG3_v1_genomic.fna \
-pe1 R1_01.fastq \
-pe2 R2_01.fastq \
-mapper minimap2 -preset sr \
-bam -rmo
samtools index minimap2.BAM/R1_01.fastq.GCA_004325075.1_ILG3_v1_genomic.fna.minimap2.bam ### Average PE depth  36.79X

samtools view -b minimap2.BAM/R1_01.fastq.GCA_004325075.1_ILG3_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I  > theta.ILG3.c4.36x.txt; mlRho -m 1000 -M 1005  >> theta.ILG3.c4.36x.txt

samtools view -b minimap2.BAM/R1_01.fastq.GCA_004325075.1_ILG3_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 8
mlRho -M 0 -I  > theta.ILG3.c8.36x.txt; mlRho -m 1000 -M 1005  >> theta.ILG3.c8.36x.txt

samtools view -b minimap2.BAM/R1_01.fastq.GCA_004325075.1_ILG3_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 16
mlRho -M 0 -I  > theta.ILG3.c16.36x.txt; mlRho -m 1000 -M 1005  >> theta.ILG3.c16.36x.txt

### on MiniMe
/media/Data_3/Karen/MlRho/

### Data from NCBI
mkdir Ordospora_colligata_GBEP
cd Ordospora_colligata_GBEP/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/324/935/GCA_004324935.1_Ordospora_GBEP_v1/GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna.gz
gunzip GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna.gz

mkdir Ordospora_colligata_FISK
cd Ordospora_colligata_FISK/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/325/055/GCA_004325055.1_Ordospora_FISK_v1/GCA_004325055.1_Ordospora_FISK_v1_genomic.fna.gz
gunzip GCA_004325055.1_Ordospora_FISK_v1_genomic.fna.gz

mkdir Ordospora_colligata_NOV7
cd Ordospora_colligata_NOV7
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/324/945/GCA_004324945.1_Ordospora_NOV7_v1/GCA_004324945.1_Ordospora_NOV7_v1_genomic.fna.gz
gunzip GCA_004324945.1_Ordospora_NOV7_v1_genomic.fna.gz

mkdir Ordospora_colligata_OC4
cd Ordospora_colligata_OC4/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/803/265/GCF_000803265.1_ASM80326v1/GCF_000803265.1_ASM80326v1_genomic.fna.gz
gunzip GCF_000803265.1_ASM80326v1_genomic.fna.gz

##### GBEP #####
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
GB-EP-1_R1.fastq.gz GB-EP-1_R2.fastq.gz \
GBEP1_R1.PE.trimmo.fastq GBEP1_R1.SE.trimmo.fastq \
GBEP1_R2.PE.trimmo.fastq GBEP1_R2.SE.trimmo.fastq \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:100

### Read mapping in PE mode with VarScan2 to check for diploid SNPs
get_SNPs.pl -fa GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna \
-pe1 GBEP1_R1.PE.trimmo.fastq \
-pe2 GBEP1_R2.PE.trimmo.fastq \
-mapper minimap2 -preset sr \
-bam \
-var /opt/varscan/VarScan.v2.4.3.jar \
-mvf 0.2

samtools index minimap2.BAM/GBEP1_R1.PE.trimmo.fastq.GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna.minimap2.bam ### Average PE depth 47.35X

samtools view -b minimap2.BAM/GBEP1_R1.PE.trimmo.fastq.GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I > theta.GBEP.c4.PE.txt; mlRho -m 1000 -M 1005 >> theta.GBEP.c4.PE.txt
# WARNING: Lower confidence limit of \Delta cannot be estimated; setting it to -1.
# Try again in SE mode to check if the issue can be resolved

### Read mapping in SE mode
cat GBEP1_R1.PE.trimmo.fastq GBEP1_R1.SE.trimmo.fastq GBEP1_R2.PE.trimmo.fastq GBEP1_R2.SE.trimmo.fastq > GBEP1.trimmo.fastq

get_SNPs.pl -fa GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna \
-fq GBEP1.trimmo.fastq \
-mapper minimap2 -preset sr \
-rmo -bam
samtools index minimap2.BAM/GBEP1.trimmo.fastq.GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna.minimap2.bam ### Average SE depth 56.96X

samtools view -b minimap2.BAM/GBEP1.trimmo.fastq.GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I > theta.GBEP.c4.SE.txt; mlRho -m 1000 -M 1005 >> theta.GBEP.c4.SE.txt
# WARNING: Lower confidence limit of \Delta cannot be estimated; setting it to -1.
# Same warning; let's try with another dataset
# Looks like not enough heterozigosity;
# compatible with Selman et al. Eukaryot Cell. 2013 Apr;12(4):496-502. doi: 10.1128/EC.00307-12.

samtools view -b minimap2.BAM/GBEP1.trimmo.fastq.GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 8
mlRho -M 0 -I > theta.GBEP.c8.SE.txt; mlRho -m 1000 -M 1005 >> theta.GBEP.c8.SE.txt

samtools view -b minimap2.BAM/GBEP1.trimmo.fastq.GCA_004324935.1_Ordospora_GBEP_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 16
mlRho -M 0 -I > theta.GBEP.c16.SE.txt; mlRho -m 1000 -M 1005 >> theta.GBEP.c16.SE.txt

##### FISK #####
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
FI-SK-17-1.R1.fastq.gz FI-SK-17-1.R2.fastq.gz \
FISK_R1.PE.trimmo.fastq FISK_R1.SE.trimmo.fastq \
FISK_R2.PE.trimmo.fastq FISK_R2.SE.trimmo.fastq \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:100

### Read mapping in SE mode
cat FISK_R1.PE.trimmo.fastq FISK_R1.SE.trimmo.fastq FISK_R2.PE.trimmo.fastq FISK_R2.SE.trimmo.fastq > FISK.trimmo.fastq

get_SNPs.pl -fa GCA_004325055.1_Ordospora_FISK_v1_genomic.fna \
-fq FISK.trimmo.fastq \
-mapper minimap2 -preset sr \
-rmo -bam
samtools index minimap2.BAM/FISK.trimmo.fastq.GCA_004325055.1_Ordospora_FISK_v1_genomic.fna.minimap2.bam ### Average SE depth 8.01X

samtools view -b minimap2.BAM/FISK.trimmo.fastq.GCA_004325055.1_Ordospora_FISK_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I > theta.FISK.c4.SE.txt; mlRho -m 1000 -M 1005 >> theta.FISK.c4.SE.txt
# WARNING: Lower confidence limit of \Delta cannot be estimated; setting it to -1.
# Same warning; let's try with another dataset

samtools view -b minimap2.BAM/FISK.trimmo.fastq.GCA_004325055.1_Ordospora_FISK_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 8
mlRho -M 0 -I > theta.FISK.c8.SE.txt; mlRho -m 1000 -M 1005 >> theta.FISK.c8.SE.txt

samtools view -b minimap2.BAM/FISK.trimmo.fastq.GCA_004325055.1_Ordospora_FISK_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 16
mlRho -M 0 -I > theta.FISK.c16.SE.txt; mlRho -m 1000 -M 1005 >> theta.FISK.c16.SE.txt

##### NOV7 #####
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
NO-V-7_R1.fastq.gz NO-V-7_R2.fastq.gz \
NOV7_R1.PE.trimmo.fastq NOV7_R1.SE.trimmo.fastq \
NOV7_R2.PE.trimmo.fastq NOV7_R2.SE.trimmo.fastq \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:100

### Read mapping in SE mode
cat NOV7_R1.PE.trimmo.fastq NOV7_R1.SE.trimmo.fastq NOV7_R2.PE.trimmo.fastq NOV7_R2.SE.trimmo.fastq > NOV7.trimmo.fastq

get_SNPs.pl -fa GCA_004324945.1_Ordospora_NOV7_v1_genomic.fna \
-fq NOV7.trimmo.fastq \
-mapper minimap2 -preset sr \
-rmo -bam
samtools index minimap2.BAM/NOV7.trimmo.fastq.GCA_004324945.1_Ordospora_NOV7_v1_genomic.fna.minimap2.bam ### Average SE depth 54.25X

samtools view -b minimap2.BAM/NOV7.trimmo.fastq.GCA_004324945.1_Ordospora_NOV7_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I > theta.NOV7.c4.SE.txt; mlRho -m 1000 -M 1005 >> theta.NOV7.c4.SE.txt

samtools view -b minimap2.BAM/NOV7.trimmo.fastq.GCA_004324945.1_Ordospora_NOV7_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 8
mlRho -M 0 -I > theta.NOV7.c8.SE.txt; mlRho -m 1000 -M 1005 >> theta.NOV7.c8.SE.txt

samtools view -b minimap2.BAM/NOV7.trimmo.fastq.GCA_004324945.1_Ordospora_NOV7_v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 16
mlRho -M 0 -I > theta.NOV7.c16.SE.txt; mlRho -m 1000 -M 1005 >> theta.NOV7.c16.SE.txt

##### OC4 ##### -> Data from Pombert et al. mBio. 2015, read length 101 bp
java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
110909_SN132_A_L004_R1_GSR-3b.fastq.gz 110909_SN132_A_L004_R2_GSR-3b.fastq.gz \
OC4_R1.PE.trimmo.fastq OC4_R1.SE.trimmo.fastq \
OC4_R2.PE.trimmo.fastq OC4_R2.SE.trimmo.fastq \
ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:75 ## lowering minlen since reads are 101 bp

cat OC4_R1.PE.trimmo.fastq OC4_R1.SE.trimmo.fastq OC4_R2.PE.trimmo.fastq OC4_R2.SE.trimmo.fastq > OC4.trimmo.fastq

get_SNPs.pl -fa GCF_000803265.1_ASM80326v1_genomic.fna \
-fq OC4.trimmo.fastq \
-mapper minimap2 -preset sr \
-rmo -bam
samtools index minimap2.BAM/OC4.trimmo.fastq.GCF_000803265.1_ASM80326v1_genomic.fna.minimap2.bam ### Average SE depth 420.57X

samtools view -b minimap2.BAM/OC4.trimmo.fastq.GCF_000803265.1_ASM80326v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 4
mlRho -M 0 -I > theta.OC4.c4.SE.txt; mlRho -m 1000 -M 1005 >> theta.OC4.c4.SE.txt

samtools view -b minimap2.BAM/OC4.trimmo.fastq.GCF_000803265.1_ASM80326v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 8
mlRho -M 0 -I > theta.OC4.c8.SE.txt; mlRho -m 1000 -M 1005 >> theta.OC4.c8.SE.txt

samtools view -b minimap2.BAM/OC4.trimmo.fastq.GCF_000803265.1_ASM80326v1_genomic.fna.minimap2.bam |
samtools mpileup - |
cut -f 2,5 |
awk -f /opt/MlRho_2.9/bam2pro.awk |
formatPro -c 16
mlRho -M 0 -I > theta.OC4.c16.SE.txt; mlRho -m 1000 -M 1005 >> theta.OC4.c16.SE.txt
