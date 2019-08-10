#!/usr/bin/bash

echo "Working on STREPTOPHYCEAE..."
## KLEBS vs WORLD
mkdir KLEBSvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying KLEBS against $DB...";\
blastp -num_threads 10 -query FASTA/Klebsormidium_v1_protein.faa -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out KLEBSvs/KLEBSvs$DB.blastp.6; \
done

echo "Working on PRASINOPHYCEAE..."
## CCMP vs WORLD
mkdir CCMPvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying CCMP against $DB...";\
blastp -num_threads 10 -query FASTA/CCMP1205_proteins.fasta -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out CCMPvs/CCMPvs$DB.blastp.6; \
done

## BATHY vs WORLD
mkdir BATHYvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying BATHY against $DB...";\
blastp -num_threads 10 -query FASTA/bathy_PROT.fasta -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out BATHYvs/BATHYvs$DB.blastp.6; \
done

## MPUSI vs WORLD
mkdir MPUSIvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying MPUSI against $DB...";\
blastp -num_threads 10 -query FASTA/MicpuC3v2_GeneCatalog_proteins_20160125.aa.fasta.newheaders -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out MPUSIvs/MPUSIvs$DB.blastp.6; \
done

## MCOMMODA vs WORLD
mkdir MCOMMODAvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying MCOMMODA against $DB...";\
blastp -num_threads 10 -query FASTA/mcommodaRCC299.fasta -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out MCOMMODAvs/MCOMMODAvs$DB.blastp.6; \
done

## OLUCI vs WORLD
mkdir OLUCIvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying OLUCI against $DB...";\
blastp -num_threads 10 -query FASTA/O.lucimarinus.FM.aa.fasta.newheaders -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out OLUCIvs/OLUCIvs$DB.blastp.6; \
done

## OTAURI vs WORLD
mkdir OTAURIvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying OTAURI against $DB...";\
blastp -num_threads 10 -query FASTA/Otauri_version_050606_protein.faa -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out OTAURIvs/OTAURIvs$DB.blastp.6; \
done

echo "Working on TREBOUXIOPHYCEAE..."
## AUXENO vs WORLD
mkdir AUXENOvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying AUXENO against $DB...";\
blastp -num_threads 10 -query FASTA/Auxenochlorella_v1_protein.faa -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out AUXENOvs/AUXENOvs$DB.blastp.6; \
done

## CHLORELLA vs WORLD
mkdir CHLORELLAvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying CHLORELLA against $DB...";\
blastp -num_threads 10 -query FASTA/Chlorella_1.0_protein.faa -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out CHLORELLAvs/CHLORELLAvs$DB.blastp.6; \
done

## COCCO vs WORLD
mkdir COCCOvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying COCCO against $DB...";\
blastp -num_threads 10 -query FASTA/Coccomyxa_v2.0_protein.faa -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out COCCOvs/COCCOvs$DB.blastp.6; \
done

## HELICO vs WORLD
mkdir HELICOvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying HELICO against $DB";\
blastp -num_threads 10 -query FASTA/Helico_v1.0_protein.faa -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out HELICOvs/HELICOvs$DB.blastp.6; \
done

echo "Working on CHLOROPHYCEAE..."
## CHLAMY vs WORLD
mkdir CHLAMYvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying CHLAMY against $DB...";\
blastp -num_threads 10 -query FASTA/Creinhardtii_281_v5.5.protein.fa -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out CHLAMYvs/CHLAMYvs$DB.blastp.6; \
done

## GONIUM vs WORLD
mkdir GONIUMvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying GONIUM against $DB...";\
blastp -num_threads 10 -query FASTA/Gonium_v1_protein.faa -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out GONIUMvs/GONIUMvs$DB.blastp.6; \
done

## VOLVOX vs WORLD
mkdir VOLVOXvs
for DB in {AUXENO,BATHY,CCMP,CHLAMY,CHLORELLA,COCCO,GONIUM,HELICO,KLEBS,MCOMMODA,MPUSI,OLUCI,OTAURI,VOLVOX}; do date; echo "Querying VOLVOX against $DB...";\
blastp -num_threads 10 -query FASTA/volvox.faa -db DB/$DB -evalue 1e-10 -culling_limit 1 -outfmt 6 -out VOLVOXvs/VOLVOXvs$DB.blastp.6; \
done
