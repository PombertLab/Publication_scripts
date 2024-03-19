#!/usr/bin/env python3

## Used for Table S4

name = 'SupplementaryMasterTable.py'
version = '0.2.0'
updated = '2022-03-20'

from argparse import ArgumentParser
from re import findall
from os.path import basename

GetOptions = ArgumentParser()
GetOptions.add_argument("-l","--loci")
GetOptions.add_argument("-s","--swissprot")
GetOptions.add_argument("-i","--interpro")
GetOptions.add_argument("-t","--trembl")
GetOptions.add_argument("-q","--quego",nargs='+')
GetOptions.add_argument("-p","--pombase",nargs='+')
GetOptions.add_argument("-g","--gesamt")
GetOptions.add_argument("-m","--mican")

args = GetOptions.parse_args()
loci_file = args.loci
interpro_file = args.interpro
swissprot_file = args.swissprot
trembl_file = args.trembl
quego_files = args.quego
pombase_files = args.pombase
gesamt_file = args.gesamt
mican_file = args.mican

max_hits = 5

############################################################
## Get loci for table
############################################################

loci = {}
annotations = {}
ordered = []
prefix = ""
LOCI = open(loci_file,'r')
for line in LOCI:
	line = line.strip()
	if line != '':
		if line[0:3] != '###':
			if line[0] == '#':
				prefix = line[3:]
			else:
				annot,locus = line.split("\t")
				locus = f"{prefix}{locus}"
				loci[locus] = {
					'annotation':annot,
					'homology':{
						'Seq':{
							'InterProScan':{
								'CDD':{
									0:{
										'score':"---",
										'ref':"---",
										'type':"E-value",
										'desc':"---"
									}
								},
								'PANTHER':{
									0:{
										'score':"---",
										'ref':"---",
										'type':"E-value",
										'desc':"---"
									}
								},
								'Pfam':{
									0:{
										'score':"---",
										'ref':"---",
										'type':"E-value",
										'desc':"---"
									}
								}
							},
							'DIAMOND':{
								'SwissProt':{
									0:{
										'score':"---",
										'ref':"---",
										'type':"E-value",
										'desc':"---"
									}
								},
								'trEMBL':{
									0:{
										'score':"---",
										'ref':"---",
										'type':"E-value",
										'desc':"---"
									}
								}
							}
						},
						'Str':{
							'QueGO':{},
							'FoldSeek':{},
							'GESAMT (MICAN)':{}
						}
					}
				}
				if annot not in annotations.keys():
					annotations[annot] = []
					ordered.append(annot)
				annotations[annot].append(locus)
LOCI.close()

############################################################
## Get SwissProt results
############################################################

SWISS = open(swissprot_file,'r')
for line in SWISS:
	if line[0] != "#":
		locus,evalue,ref,desc= line.strip().split("\t")
		if evalue != 'N/A':
			if locus in loci.keys():
				loci[locus]['homology']['Seq']['DIAMOND']['SwissProt'] = {
					1:{
						'score':evalue,
						'ref':ref,
						'desc':desc,
						'type':'E-value'
					}
				}
SWISS.close()

# ############################################################
# ## Get trEMBL results
# ############################################################

TREMBL = open(trembl_file,'r')
for line in TREMBL:
	if line[0] != "#":
		locus,evalue,ref,desc= line.strip().split("\t")
		if evalue != 'N/A':
			if locus in loci.keys():
				loci[locus]['homology']['Seq']['DIAMOND']['trEMBL'] = {
					1:{
						'score':evalue,
						'ref':ref,
						'desc':desc,
						'type':'E-value'
					}
				}
TREMBL.close()

############################################################
## Get InterPro results
############################################################

sequences = {}

desired_dbs = ['PANTHER','Pfam','CDD']
INTER = open(interpro_file,'r')
db_count = 1
for line in INTER:
	
	line = line.strip()
	data = line.split("\t")
	locus = data[0]
	
	if locus in loci.keys():
		
		db = data[3]

		if db in desired_dbs:

			if locus not in sequences.keys():
				sequences[locus] = {}

			if db_count in loci[locus]['homology']['Seq']['InterProScan'][db].keys():
				db_count += 1
			else:
				db_count = 1

			ref = data[4]
			desc = data[5]
			score = data[8]

			if ref not in sequences[locus].keys():

				if db_count < max_hits:

					loci[locus]['homology']['Seq']['InterProScan'][db][db_count] = {
						'score':score,
						'ref':ref,
						'desc':desc,
						'type':'E-value',
					}

					sequences[locus][ref] = True

INTER.close()

############################################################
## Get QueGO results
############################################################

structures = {}
sequences = {}

for quego_file in quego_files:

	QUEGO = open(quego_file,'r')
	go = findall("QueGO_(\d+)",basename(quego_file).split(".")[0])[0]
	annotation = ""
	
	seq_count = 1
	str_count = 1

	for locus in loci.keys():
		loci[locus]['homology']['Seq']['DIAMOND'][f'UniProt:GO{go}'] = {}
		loci[locus]['homology']['Str']['QueGO'][f'UniProt:GO{go}'] = {}

	for line in QUEGO:
		
		line = line.strip()
		
		if line != '':
			
			if line[0:3] != "###":
				
				if line[0:2] == "##":
					
					annotation = line.split("\t")[0][3:]
				
				else:
					
					locus,annot,ss,sr,fs,fr,fd,gs,gr,gd = line.split("\t")[0:10]

					if locus in loci.keys():
					
						if locus not in structures.keys():
							structures[locus] = {}
							sequences[locus] = {}

						if fs != '-' and fr not in structures[locus].keys():

							if str_count in loci[locus]['homology']['Str']['QueGO'][f'UniProt:GO{go}'].keys():
								str_count += 1
							else:
								str_count = 1

							if str_count < max_hits:
							
								loci[locus]['homology']['Str']['QueGO'][f'UniProt:GO{go}'][str_count] = {
									'score':fs,
									'ref':fr,
									'type':'TM-Score',
									'desc':annotation
								}

								structures[locus][fr] = True

						if ss != '-'and sr not in sequences[locus].keys():
							
							if seq_count in loci[locus]['homology']['Seq']['DIAMOND'][f'UniProt:GO{go}'].keys():
								seq_count += 1
							else:
								seq_count = 1

							if seq_count < max_hits:

								loci[locus]['homology']['Seq']['DIAMOND'][f'UniProt:GO{go}'][seq_count] = {
									'score':ss,
									'ref':sr,
									'type':'E-value',
									'desc':annotation
								}

								sequences[locus][sr] = True
	QUEGO.close()

############################################################
## Get PomBase results
############################################################

structures = {}
sequences = {}

for pombase_file in pombase_files:
	POMBASE = open(pombase_file,'r')
	seq_count = 1
	str_count = 1
	go = ""
	for line in POMBASE:
		line = line.strip()
		if line != '':
			if line[0:5] == '#####':
				go = findall("\(GO:(\d+)\)",line)[0]

				for locus in loci.keys():
					loci[locus]['homology']['Seq']['DIAMOND'][f'PomBase:GO{go}'] = {
					}
					loci[locus]['homology']['Str']['FoldSeek'][f'PomBase:GO{go}'] = {
					}

			else:
				if line[0] != '#':

					gene,desc,locus,annot,be,fe,tm = line.split("\t")

					if locus in loci.keys():

						if locus not in structures.keys():
							structures[locus] = {}
							sequences[locus] = {}

						if tm != '---' and gene not in structures[locus].keys():

							if str_count in loci[locus]['homology']['Str']['FoldSeek'][f'PomBase:GO{go}'].keys():
								str_count += 1
							else:
								str_count = 1
							
							if str_count < max_hits:

								loci[locus]['homology']['Str']['FoldSeek'][f'PomBase:GO{go}'][str_count] = {
									'score':tm,
									'ref':gene,
									'type':'TM-score',
									'desc':desc
								}

								structures[locus][gene] = True

						if be != '---' and gene not in sequences[locus].keys():

							if seq_count in loci[locus]['homology']['Seq']['DIAMOND'][f'PomBase:GO{go}'].keys():
								seq_count += 1
							else:
								seq_count = 1
							
							if seq_count < max_hits:

								loci[locus]['homology']['Seq']['DIAMOND'][f'PomBase:GO{go}'][seq_count] = {
									'score':be,
									'ref':gene,
									'type':'E-value',
									'desc':desc
								}

								sequences[locus][gene] = True
	POMBASE.close()

############################################################
## Get 3DFI results
############################################################

GESAMT = open(gesamt_file,'r')
str_count = 1
gesamt_pdb_links = {}

for locus in loci.keys():
	
	loci[locus]['homology']['Str']['GESAMT (MICAN)']['RCSB'] = {}

	gesamt_pdb_links[locus] = {}

for line in GESAMT:

	line = line.strip()


	if line != '':

		if line[0] != '#':

			mod,pred,pdb,chain,qscore,stat1,rmsd,start,end,filename,desc = line.split("\t")
			locus = findall("^(\w+)-m\d+",mod)[0]

			if locus in loci.keys():

				if str_count in loci[locus]['homology']['Str']['GESAMT (MICAN)']['RCSB'].keys():
					str_count += 1
				else:
					str_count = 1
				
				if str_count < max_hits:

					ref = f"{pdb}_{chain}"

					loci[locus]['homology']['Str']['GESAMT (MICAN)']['RCSB'][str_count] = {
						'score':qscore,
						'ref':ref,
						'type':'Q-score (TM-score)',
						'desc':desc
					}

					gesamt_pdb_links[locus][ref] = str_count


GESAMT.close()

############################################################
## Write output
############################################################
count = 0
MICAN = open(mican_file,'r')
for line in MICAN:
	line = line.strip()
	if line != '':
		if line[0] != '#':
			locus,pred,pdb,chain,stm,tm = line.split("\t")[0:6]
			ref = f"{pdb}_{chain}"
			locus = findall("^(\w+)-m\d+",locus)[0]
			if ref in gesamt_pdb_links[locus].keys():
				str_count = gesamt_pdb_links[locus][f"{pdb}_{chain}"]
				qscore = loci[locus]['homology']['Str']['GESAMT (MICAN)']['RCSB'][str_count]['score']
				ref = loci[locus]['homology']['Str']['GESAMT (MICAN)']['RCSB'][str_count]['ref']
				loci[locus]['homology']['Str']['GESAMT (MICAN)']['RCSB'][str_count]['score'] = f"{qscore} ({tm})"
MICAN.close()

############################################################
## Write output
############################################################

## Final output file should look like the following:
## Locus Tag - Annotation
## HomType	Method	Score	Database	Reference	Description

OUT = open("master_table.tsv",'w')
OUT.write(f"## Locus\tHomologyTool\tDatabase\tHomologyScore\tHomologyScoreType\tHomologyHit\tDescription\n")
for annot in ordered:
	# print(annot)
	OUT.write(f"### {annot}\n")
	for locus in annotations[annot]:
		# print(f"\t{locus}")
		# print(loci[locus])
		for hom_type in loci[locus]['homology']:
			# print(f"\t\t{hom_type}")
			OUT.write("\n")
			for tool in loci[locus]['homology'][hom_type]:
				# print(f"\t\t\t{tool}")
				for db in loci[locus]['homology'][hom_type][tool].keys():
					# print(f"\t\t\t\t{db}")
					if hom_type == 'Seq':
						for result in sorted(loci[locus]['homology'][hom_type][tool][db],key=lambda x: loci[locus]['homology'][hom_type][tool][db][x]['score']):
							# print(f"\t\t\t\t\t{result}")
							# print(f"\t\t\t\t\t\t{loci[locus]['homology'][hom_type][tool][db][result]}")
							if len(loci[locus]['homology'][hom_type][tool][db].keys()) > 1 and result == 0:
								continue
							score = loci[locus]['homology'][hom_type][tool][db][result]['score']
							score_type = loci[locus]['homology'][hom_type][tool][db][result]['type']
							ref = loci[locus]['homology'][hom_type][tool][db][result]['ref']
							desc = loci[locus]['homology'][hom_type][tool][db][result]['desc']
							OUT.write(f"{locus}\t{tool}\t{db}\t{score}\t{score_type}\t{ref}\t{desc}\n")
					elif hom_type == 'Str':
						for result in sorted(loci[locus]['homology'][hom_type][tool][db],key=lambda x: loci[locus]['homology'][hom_type][tool][db][x]['score'],reverse=True):
							# print(f"\t\t\t\t\t{result}")
							# print(f"\t\t\t\t\t\t{loci[locus]['homology'][hom_type][tool][db][result]}")
							if len(loci[locus]['homology'][hom_type][tool][db].keys()) > 1 and result == 0:
								continue
							score = loci[locus]['homology'][hom_type][tool][db][result]['score']
							score_type = loci[locus]['homology'][hom_type][tool][db][result]['type']
							ref = loci[locus]['homology'][hom_type][tool][db][result]['ref']
							desc = loci[locus]['homology'][hom_type][tool][db][result]['desc']
							OUT.write(f"{locus}\t{tool}\t{db}\t{score}\t{score_type}\t{ref}\t{desc}\n")
		OUT.write("\n")
OUT.close()
