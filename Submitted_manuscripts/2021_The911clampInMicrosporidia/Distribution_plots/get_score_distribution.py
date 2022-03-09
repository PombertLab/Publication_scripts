#!/usr/bin/python
## Pombert Lab 2022

name = 'get_score_distribution.py'
version = '0.4'
updated = '2022-03-09'

import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse as ap
from sys import argv

usage = f"""
NAME		{name}
VERSION		{version}
UPDATED		{updated}
SYNOPSIS	This script takes the master table produced by create_master_table.pl and plots various
		distributions or comparisons.

USAGE		{name} \\
			  -plddt \\
			  -m MASTER_TABLE.tsv

OPTIONS
-m (--master_table)	Master table produced by create_master_table.pl
--plddt		Plot the pLDDT distribution
--q		Plot the Q-score distribution
--tm	Plot the TM-score distribution
--voro	Plot the voroCNN distribution
--tmvq		Plot the TM- vs Q-scores to visualize the relationship between the two
--pdvvo	Plot the pLDDT vs voroCNN score
"""

if(len(argv) < 2):
	exit(f"\n{usage}\n")

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}

mpl.rcParams.update(new_rc_params)

args = ap.ArgumentParser()

args.add_argument("-m","--master_table",required=True)
args.add_argument("--plddt",action="store_true")
args.add_argument("--q",action="store_true")
args.add_argument("--tm",action="store_true")
args.add_argument("--voro",action="store_true")
args.add_argument("--tmvq",action="store_true")
args.add_argument("--pdvvo",action="store_true")

master_file = args.master_table
plddt = args.plddt
q = args.q
tm = args.tm
voro = args.voro
tmvq = args.tmvq
pdvvo = args.pdvvo

if(tm):

	## All proteins

	tm_score = []
	a_tmscore = []
	r_tmscore = []

	all_proteins = 0

	a_tmscore_greater_than_90 = 0
	a_tmscore_greater_than_70 = 0
	a_tmscore_greater_than_50 = 0
	a_tmscore_less_than_50 = 0

	r_tmscore_greater_than_90 = 0
	r_tmscore_greater_than_70 = 0
	r_tmscore_greater_than_50 = 0
	r_tmscore_less_than_50 = 0

	## Known proteins
	k_a_tmscore = []
	k_r_tmscore = []

	known_proteins = 0
	known_proteins_w_alpha = 0
	known_proteins_w_raptor = 0

	k_a_tmscore_greater_than_90 = 0
	k_a_tmscore_greater_than_70 = 0
	k_a_tmscore_greater_than_50 = 0
	k_a_tmscore_less_than_50 = 0

	k_r_tmscore_greater_than_90 = 0
	k_r_tmscore_greater_than_70 = 0
	k_r_tmscore_greater_than_50 = 0
	k_r_tmscore_less_than_50 = 0

	## Hypothetical proteins
	h_a_tmscore = []
	h_r_tmscore = []

	hypo_proteins = 0
	hypo_proteins_w_alpha = 0
	hypo_proteins_w_raptor = 0

	h_a_tmscore_greater_than_90 = 0
	h_a_tmscore_greater_than_70 = 0
	h_a_tmscore_greater_than_50 = 0
	h_a_tmscore_less_than_50 = 0

	h_r_tmscore_greater_than_90 = 0
	h_r_tmscore_greater_than_70 = 0
	h_r_tmscore_greater_than_50 = 0
	h_r_tmscore_less_than_50 = 0

	with open (argv[2],"r") as FILE:
		for line in FILE:

			line = line.replace("\n",'')

			## Parse data for non-headers and divider lines
			if line != "" and line[0:5] != "Locus":
				data = line.split("\t")
				a_TM_score = data[8]
				r_TM_score = data[17]

				## Get data only for the best hits
				if(data[0] != ''):

					all_proteins += 1
					
					## Keep track of how many proteins are known vs unknown
					if(data[1] == "hypothetical protein"):
						hypo_proteins += 1
					elif("Uncharacterized" in data[1]):
						hypo_proteins += 1
					elif("HYPOTHETICAL" in data[1]):
						hypo_proteins += 1
					elif("UPF" in data[1]):
						hypo_proteins += 1
					else:
						known_proteins += 1
						

					## If there is a TM score for AlphaFold, store it, and count what quality it falls under
					if(a_TM_score != '---' and a_TM_score != ''):

						a_TM_score = float(data[8])
						a_tmscore.append(a_TM_score)
						
						if(a_TM_score >= .90):
							a_tmscore_greater_than_90 += 1
							a_tmscore_greater_than_70 += 1
							a_tmscore_greater_than_50 += 1
						elif(a_TM_score >= .70):
							a_tmscore_greater_than_70 += 1
							a_tmscore_greater_than_50 += 1
						elif(a_TM_score >= .50):
							a_tmscore_greater_than_50 += 1
						elif(a_TM_score < .50):
							a_tmscore_less_than_50 += 1
						
						if(data[1] == "hypothetical protein" or "Uncharacterized" in data[1] or "HYPOTHETICAL" in data[1] or "UPF" in data[1]):
							h_a_tmscore.append(a_TM_score)
							hypo_proteins_w_alpha += 1
							if(a_TM_score >= .90):
								h_a_tmscore_greater_than_90 += 1
								h_a_tmscore_greater_than_70 += 1
								h_a_tmscore_greater_than_50 += 1
							elif(a_TM_score >= .70):
								h_a_tmscore_greater_than_70 += 1
								h_a_tmscore_greater_than_50 += 1
							elif(a_TM_score >= .50):
								h_a_tmscore_greater_than_50 += 1
							elif(a_TM_score < .50):
								h_a_tmscore_less_than_50 += 1
						else:
							k_a_tmscore.append(a_TM_score)
							known_proteins_w_alpha += 1
							if(a_TM_score >= .90):
								k_a_tmscore_greater_than_90 += 1
								k_a_tmscore_greater_than_70 += 1
								k_a_tmscore_greater_than_50 += 1
							elif(a_TM_score >= .70):
								k_a_tmscore_greater_than_70 += 1
								k_a_tmscore_greater_than_50 += 1
							elif(a_TM_score >= .50):
								k_a_tmscore_greater_than_50 += 1
							elif(a_TM_score < .50):
								k_a_tmscore_less_than_50 += 1

					## If there is a TM score for AlphaFold, store it, and count what quality it falls under
					if(r_TM_score != '---' and r_TM_score != ''):

						r_TM_score = float(data[17])
						r_tmscore.append(r_TM_score)
						
						if(r_TM_score >= .90):
							r_tmscore_greater_than_90 += 1
							r_tmscore_greater_than_70 += 1
							r_tmscore_greater_than_50 += 1
						elif(r_TM_score >= .70):
							r_tmscore_greater_than_70 += 1
							r_tmscore_greater_than_50 += 1
						elif(r_TM_score >= .50):
							r_tmscore_greater_than_50 += 1
						elif(r_TM_score < .50):
							r_tmscore_less_than_50 += 1
						
						if(data[1] == "hypothetical protein" or "Uncharacterized" in data[1] or "HYPOTHETICAL" in data[1] or "UPF" in data[1]):
							h_r_tmscore.append(r_TM_score)
							hypo_proteins_w_raptor += 1
							if(r_TM_score >= .90):
								h_r_tmscore_greater_than_90 += 1
								h_r_tmscore_greater_than_70 += 1
								h_r_tmscore_greater_than_50 += 1
							elif(r_TM_score >= .70):
								h_r_tmscore_greater_than_70 += 1
								h_r_tmscore_greater_than_50 += 1
							elif(r_TM_score >= .50):
								h_r_tmscore_greater_than_50 += 1
							elif(r_TM_score < .50):
								h_r_tmscore_less_than_50 += 1
						else:
							k_r_tmscore.append(r_TM_score)
							known_proteins_w_raptor += 1
							if(r_TM_score >= .90):
								k_r_tmscore_greater_than_90 += 1
								k_r_tmscore_greater_than_70 += 1
								k_r_tmscore_greater_than_50 += 1
							elif(r_TM_score >= .70):
								k_r_tmscore_greater_than_70 += 1
								k_r_tmscore_greater_than_50 += 1
							elif(r_TM_score >= .50):
								k_r_tmscore_greater_than_50 += 1
							elif(r_TM_score < .50):
								k_r_tmscore_less_than_50 += 1

	a_header = f'''
	TM score distribution for AlphaFold Predictions
	'''

	a_header = a_header.replace("\t","")

	r_header = f'''
	TM score distribution for RaptorX Predictions
	'''

	# Plot text for hypothetical proteins
	h_atext = fr'''
	3D Homologs Found: {hypo_proteins_w_alpha}
	$\geq$0.90   {'{:4d}'.format(h_a_tmscore_greater_than_90)} ({'{:5.2f}'.format(h_a_tmscore_greater_than_90/len(h_a_tmscore)*100)}%)
	$\geq$0.70   {'{:4d}'.format(h_a_tmscore_greater_than_70)} ({'{:5.2f}'.format(h_a_tmscore_greater_than_70/len(h_a_tmscore)*100)}%)
	$\geq$0.50   {'{:4d}'.format(h_a_tmscore_greater_than_50)} ({'{:5.2f}'.format(h_a_tmscore_greater_than_50/len(h_a_tmscore)*100)}%)
	$<$0.50   {'{:4d}'.format(h_a_tmscore_less_than_50)} ({'{:5.2f}'.format((h_a_tmscore_less_than_50/len(h_a_tmscore)*100))}%)
	'''

	h_atext = h_atext.replace("\t","")

	h_rtext = fr'''
	3D Homologs Found: {hypo_proteins_w_raptor}
	$\geq$0.90   {'{:4d}'.format(h_r_tmscore_greater_than_90)} ({'{:5.2f}'.format(h_r_tmscore_greater_than_90/len(h_r_tmscore)*100)}%)
	$\geq$0.70   {'{:4d}'.format(h_r_tmscore_greater_than_70)} ({'{:5.2f}'.format(h_r_tmscore_greater_than_70/len(h_r_tmscore)*100)}%)
	$\geq$0.50   {'{:4d}'.format(h_r_tmscore_greater_than_50)} ({'{:5.2f}'.format(h_r_tmscore_greater_than_50/len(h_r_tmscore)*100)}%)
	$<$0.50   {'{:4d}'.format(h_r_tmscore_less_than_50)} ({'{:5.2f}'.format((h_r_tmscore_less_than_50/len(h_r_tmscore)*100))}%)
	'''

	h_rtext = h_rtext.replace("\t","")

	## Plot text for known proteins
	k_atext = fr'''
	3D Homologs Found: {known_proteins_w_alpha}
	$\geq$0.90   {'{:4d}'.format(k_a_tmscore_greater_than_90)} ({'{:5.2f}'.format(k_a_tmscore_greater_than_90/len(k_a_tmscore)*100)}%)
	$\geq$0.70   {'{:4d}'.format(k_a_tmscore_greater_than_70)} ({'{:5.2f}'.format(k_a_tmscore_greater_than_70/len(k_a_tmscore)*100)}%)
	$\geq$0.50   {'{:4d}'.format(k_a_tmscore_greater_than_50)} ({'{:5.2f}'.format(k_a_tmscore_greater_than_50/len(k_a_tmscore)*100)}%)
	$<$0.50   {'{:4d}'.format(k_a_tmscore_less_than_50)} ({'{:5.2f}'.format((k_a_tmscore_less_than_50/len(k_a_tmscore)*100))}%)
	'''

	k_atext = k_atext.replace("\t","")

	k_rtext = fr'''
	3D Homologs Found: {known_proteins_w_raptor}
	$\geq$0.90   {'{:4d}'.format(k_r_tmscore_greater_than_90)} ({'{:5.2f}'.format(k_r_tmscore_greater_than_90/len(k_r_tmscore)*100)}%)
	$\geq$0.70   {'{:4d}'.format(k_r_tmscore_greater_than_70)} ({'{:5.2f}'.format(k_r_tmscore_greater_than_70/len(k_r_tmscore)*100)}%)
	$\geq$0.50   {'{:4d}'.format(k_r_tmscore_greater_than_50)} ({'{:5.2f}'.format(k_r_tmscore_greater_than_50/len(k_r_tmscore)*100)}%)
	$<$0.50   {'{:4d}'.format(k_r_tmscore_less_than_50)} ({'{:5.2f}'.format((k_r_tmscore_less_than_50/len(k_r_tmscore)*100))}%)
	'''

	k_rtext = k_rtext.replace("\t","")

	# TM score distribution with hypo-protein overlay
	bins = 50
	max = 1
	fig, (alp,rpx) = plt.subplots(2,1)

	blnk = []
	all = "#6F333B"
	hypo = "#1f77b4"
	known = "#ff7f0e"

	alp.hist([h_a_tmscore,k_a_tmscore],bins=50,rwidth = .95, color=[hypo,known], range = [0,1],stacked=True, label=[f"Hypothetical Proteins: {hypo_proteins}",f"Known Proteins: {known_proteins}"])
	alp.set_xlim(0,1)
	alp.set_ylim(0,125)
	alp.set_xticks([float(max/bins)*i for i in range(bins+1)])
	alp.text(.02,110,a_header,fontsize="x-large",fontfamily="monospace")
	alp.text(.02,70,h_atext,fontsize="x-large",fontfamily="monospace",color=hypo)
	alp.text(.22,70,k_atext,fontsize="x-large",fontfamily="monospace",color=known)
	alp.legend(loc=0)

	rpx.hist([h_r_tmscore,k_r_tmscore],bins=50, rwidth = .95, color=[hypo,known], range = [0,1],stacked=True, label=[f"Hypothetical Proteins: {hypo_proteins}",f"Known Proteins: {known_proteins}"])
	rpx.set_xlim(0,1)
	rpx.set_ylim(0,125)
	rpx.set_xticks([float(max/bins)*i for i in range(bins+1)])
	rpx.text(.02,110,r_header,fontsize="x-large",fontfamily="monospace")
	rpx.text(.02,70,h_rtext,fontsize="x-large",fontfamily="monospace",color=hypo)
	rpx.text(.22,70,k_rtext,fontsize="x-large",fontfamily="monospace",color=known)
	rpx.legend(loc=0)

	plt.show()

if(plddt):
	###################################################################################################
	## pLDDT distribution graph
	###################################################################################################

	a_plddt = []

	all_proteins = 0

	a_plddt_greater_than_90 = 0
	a_plddt_greater_than_70 = 0
	a_plddt_greater_than_50 = 0
	a_plddt_less_than_50 = 0
	a_plddt_less_than_30 = 0

	## Known proteins
	k_a_plddt = []

	known_proteins = 0
	known_proteins_w_alpha = 0

	k_a_plddt_greater_than_90 = 0
	k_a_plddt_greater_than_70 = 0
	k_a_plddt_greater_than_50 = 0
	k_a_plddt_less_than_50 = 0
	k_a_plddt_less_than_30 = 0

	## Hypothetical proteins
	h_a_plddt = []

	hypo_proteins = 0
	hypo_proteins_w_alpha = 0

	h_a_plddt_greater_than_90 = 0
	h_a_plddt_greater_than_70 = 0
	h_a_plddt_greater_than_50 = 0
	h_a_plddt_less_than_50 = 0
	h_a_plddt_less_than_30 = 0

	with open (argv[2],"r") as FILE:
		for line in FILE:

			line = line.replace("\n",'')

			## Parse data for non-headers and divider lines
			if line != "" and line[0:5] != "Locus":
				data = line.split("\t")
				a_q_score = data[3]

				## Get data only for the best hits
				if(data[0] != ''):

					## If there is a q score for AlphaFold, store it, and count what quality it falls under
					if(a_q_score != '---' and a_q_score != ''):
						
						all_proteins += 1
						
						## Keep track of how many proteins are known vs unknown
						if(data[1] == "hypothetical protein"):
							hypo_proteins += 1
						elif("Uncharacterized" in data[1]):
							hypo_proteins += 1
						elif("HYPOTHETICAL" in data[1]):
							hypo_proteins += 1
						elif("UPF" in data[1]):
							hypo_proteins += 1
						else:
							known_proteins += 1

						a_q_score = float(data[3])
						a_plddt.append(a_q_score)
						
						if(a_q_score >= 90):
							a_plddt_greater_than_90 += 1
							a_plddt_greater_than_70 += 1
							a_plddt_greater_than_50 += 1
						elif(a_q_score >= 70):
							a_plddt_greater_than_70 += 1
							a_plddt_greater_than_50 += 1
						elif(a_q_score >= 50):
							a_plddt_greater_than_50 += 1
						elif(a_q_score < 50):
							a_plddt_less_than_50 += 1
						
						if(data[1] == "hypothetical protein" or "Uncharacterized" in data[1] or "HYPOTHETICAL" in data[1] or "UPF" in data[1]):
							h_a_plddt.append(a_q_score)
							hypo_proteins_w_alpha += 1
							if(a_q_score >= 90):
								h_a_plddt_greater_than_90 += 1
								h_a_plddt_greater_than_70 += 1
								h_a_plddt_greater_than_50 += 1
							elif(a_q_score >= 70):
								h_a_plddt_greater_than_70 += 1
								h_a_plddt_greater_than_50 += 1
							elif(a_q_score >= 50):
								h_a_plddt_greater_than_50 += 1
							elif(a_q_score < 50):
								h_a_plddt_less_than_50 += 1
						else:
							k_a_plddt.append(a_q_score)
							known_proteins_w_alpha += 1
							if(a_q_score >= 90):
								k_a_plddt_greater_than_90 += 1
								k_a_plddt_greater_than_70 += 1
								k_a_plddt_greater_than_50 += 1
							elif(a_q_score >= 70):
								k_a_plddt_greater_than_70 += 1
								k_a_plddt_greater_than_50 += 1
							elif(a_q_score >= 50):
								k_a_plddt_greater_than_50 += 1
							elif(a_q_score < 50):
								k_a_plddt_less_than_50 += 1
								
	a_header = f'''
	pLDDT distribution for AlphaFold Predictions
	'''

	a_header = a_header.replace("\t","")

	# Plot text for hypothetical proteins
	h_atext = fr'''
	$\geq$90   {'{:4d}'.format(h_a_plddt_greater_than_90)} ({'{:5.2f}'.format(h_a_plddt_greater_than_90/len(h_a_plddt)*100)}%)
	$\geq$70   {'{:4d}'.format(h_a_plddt_greater_than_70)} ({'{:5.2f}'.format(h_a_plddt_greater_than_70/len(h_a_plddt)*100)}%)
	$\geq$50   {'{:4d}'.format(h_a_plddt_greater_than_50)} ({'{:5.2f}'.format(h_a_plddt_greater_than_50/len(h_a_plddt)*100)}%)
	$<$50   {'{:4d}'.format(h_a_plddt_less_than_50)} ({'{:5.2f}'.format((h_a_plddt_less_than_50/len(h_a_plddt)*100))}%)
	'''

	h_atext = h_atext.replace("\t","")

	## Plot text for known proteins
	k_atext = fr'''
	$\geq$90   {'{:4d}'.format(k_a_plddt_greater_than_90)} ({'{:5.2f}'.format(k_a_plddt_greater_than_90/len(k_a_plddt)*100)}%)
	$\geq$70   {'{:4d}'.format(k_a_plddt_greater_than_70)} ({'{:5.2f}'.format(k_a_plddt_greater_than_70/len(k_a_plddt)*100)}%)
	$\geq$50   {'{:4d}'.format(k_a_plddt_greater_than_50)} ({'{:5.2f}'.format(k_a_plddt_greater_than_50/len(k_a_plddt)*100)}%)
	$<$50   {'{:4d}'.format(k_a_plddt_less_than_50)} ({'{:5.2f}'.format((k_a_plddt_less_than_50/len(k_a_plddt)*100))}%)
	'''

	k_atext = k_atext.replace("\t","")

	# q score distribution with hypo-protein overlay
	bins = 50
	max = 100
	fig, (alp) = plt.subplots(1,1)

	blnk = []
	all = "#6F333B"
	hypo = "#1f77b4"
	known = "#ff7f0e"

	alp.hist([h_a_plddt,k_a_plddt],bins=50,rwidth = .95, color=[hypo,known], range = [0,100],stacked=True, label=[f"Hypothetical Proteins: {hypo_proteins}",f"Known Proteins: {known_proteins}"])
	alp.set_xlim(0,100)
	alp.set_ylim(0,205)
	alp.set_xticks([float(max/bins)*i for i in range(0,bins+1)])
	alp.text(10,150,a_header,fontsize="x-large",fontfamily="monospace")
	alp.text(10,125,h_atext,fontsize="x-large",fontfamily="monospace",color=hypo)
	alp.text(30,125,k_atext,fontsize="x-large",fontfamily="monospace",color=known)
	alp.legend(loc=2)

	plt.show()

if(q):
	###################################################################################################
	## Q-score distribution graph
	###################################################################################################

	a_qscore = []
	r_qscore = []

	all_proteins = 0

	a_qscore_greater_than_90 = 0
	a_qscore_greater_than_70 = 0
	a_qscore_greater_than_50 = 0
	a_qscore_less_than_50 = 0
	a_qscore_less_than_30 = 0

	r_qscore_greater_than_90 = 0
	r_qscore_greater_than_70 = 0
	r_qscore_greater_than_50 = 0
	r_qscore_less_than_50 = 0
	r_qscore_less_than_30 = 0

	## Known proteins
	k_a_qscore = []
	k_r_qscore = []

	known_proteins = 0
	known_proteins_w_alpha = 0
	known_proteins_w_raptor = 0

	k_a_qscore_greater_than_90 = 0
	k_a_qscore_greater_than_70 = 0
	k_a_qscore_greater_than_50 = 0
	k_a_qscore_less_than_50 = 0
	k_a_qscore_less_than_30 = 0

	k_r_qscore_greater_than_90 = 0
	k_r_qscore_greater_than_70 = 0
	k_r_qscore_greater_than_50 = 0
	k_r_qscore_less_than_50 = 0
	k_r_qscore_less_than_30 = 0

	## Hypothetical proteins
	h_a_qscore = []
	h_r_qscore = []

	hypo_proteins = 0
	hypo_proteins_w_alpha = 0
	hypo_proteins_w_raptor = 0

	h_a_qscore_greater_than_90 = 0
	h_a_qscore_greater_than_70 = 0
	h_a_qscore_greater_than_50 = 0
	h_a_qscore_less_than_50 = 0
	h_a_qscore_less_than_30 = 0

	h_r_qscore_greater_than_90 = 0
	h_r_qscore_greater_than_70 = 0
	h_r_qscore_greater_than_50 = 0
	h_r_qscore_less_than_50 = 0
	h_r_qscore_less_than_30 = 0

	with open (argv[2],"r") as FILE:
		for line in FILE:

			line = line.replace("\n",'')

			## Parse data for non-headers and divider lines
			if line != "" and line[0:5] != "Locus":
				data = line.split("\t")
				a_q_score = data[7]
				r_q_score = data[16]

				## Get data only for the best hits
				if(data[0] != ''):

					all_proteins += 1
					
					## Keep track of how many proteins are known vs unknown
					if(data[1] == "hypothetical protein"):
						hypo_proteins += 1
					elif("Uncharacterized" in data[1]):
						hypo_proteins += 1
					elif("HYPOTHETICAL" in data[1]):
						hypo_proteins += 1
					elif("UPF" in data[1]):
						hypo_proteins += 1
					else:
						known_proteins += 1
						

					## If there is a q score for AlphaFold, store it, and count what quality it falls under
					if(a_q_score != '---' and a_q_score != ''):

						a_q_score = float(data[7])
						a_qscore.append(a_q_score)
						
						if(a_q_score >= .90):
							a_qscore_greater_than_90 += 1
							a_qscore_greater_than_70 += 1
							a_qscore_greater_than_50 += 1
							a_qscore_less_than_50 += 1
						elif(a_q_score >= .70):
							a_qscore_greater_than_70 += 1
							a_qscore_greater_than_50 += 1
							a_qscore_less_than_50 += 1
						elif(a_q_score >= .50):
							a_qscore_greater_than_50 += 1
							a_qscore_less_than_50 += 1
						elif(a_q_score >= .30):
							a_qscore_less_than_50 += 1
						elif(a_q_score < .3):
							a_qscore_less_than_30 += 1
						
						if(data[1] == "hypothetical protein" or "Uncharacterized" in data[1] or "HYPOTHETICAL" in data[1] or "UPF" in data[1]):
							h_a_qscore.append(a_q_score)
							hypo_proteins_w_alpha += 1
							if(a_q_score >= .90):
								h_a_qscore_greater_than_90 += 1
								h_a_qscore_greater_than_70 += 1
								h_a_qscore_greater_than_50 += 1
								h_a_qscore_less_than_50 += 1
							elif(a_q_score >= .70):
								h_a_qscore_greater_than_70 += 1
								h_a_qscore_greater_than_50 += 1
								h_a_qscore_less_than_50 += 1
							elif(a_q_score >= .50):
								h_a_qscore_greater_than_50 += 1
								h_a_qscore_less_than_50 += 1
							elif(a_q_score >= .30):
								h_a_qscore_less_than_50 += 1
							elif(a_q_score < .3):
								h_a_qscore_less_than_30 += 1
						else:
							k_a_qscore.append(a_q_score)
							known_proteins_w_alpha += 1
							if(a_q_score >= .90):
								k_a_qscore_greater_than_90 += 1
								k_a_qscore_greater_than_70 += 1
								k_a_qscore_greater_than_50 += 1
								k_a_qscore_less_than_50 += 1
							elif(a_q_score >= .70):
								k_a_qscore_greater_than_70 += 1
								k_a_qscore_greater_than_50 += 1
								k_a_qscore_less_than_50 += 1
							elif(a_q_score >= .50):
								k_a_qscore_greater_than_50 += 1
								k_a_qscore_less_than_50 += 1
							elif(a_q_score >= .30):
								k_a_qscore_less_than_50 += 1
							elif(a_q_score < .3):
								k_a_qscore_less_than_30 += 1

					## If there is a q score for AlphaFold, store it, and count what quality it falls under
					if(r_q_score != '---' and r_q_score != ''):

						r_q_score = float(data[16])
						r_qscore.append(r_q_score)

						if(r_q_score >= .90):
							r_qscore_greater_than_90 += 1
							r_qscore_greater_than_70 += 1
							r_qscore_greater_than_50 += 1
							r_qscore_less_than_50 += 1
						elif(r_q_score >= .70):
							r_qscore_greater_than_70 += 1
							r_qscore_greater_than_50 += 1
							r_qscore_less_than_50 += 1
						elif(r_q_score >= .50):
							r_qscore_greater_than_50 += 1
							r_qscore_less_than_50 += 1
						elif(r_q_score >= .30):
							r_qscore_less_than_50 += 1
						elif(r_q_score < .3):
							r_qscore_less_than_30 += 1
						
						if(data[1] == "hypothetical protein" or "Uncharacterized" in data[1] or "HYPOTHETICAL" in data[1] or "UPF" in data[1]):
							h_r_qscore.append(r_q_score)
							hypo_proteins_w_raptor += 1
							if(r_q_score >= .90):
								h_r_qscore_greater_than_90 += 1
								h_r_qscore_greater_than_70 += 1
								h_r_qscore_greater_than_50 += 1
								h_r_qscore_less_than_50 += 1
							elif(r_q_score >= .70):
								h_r_qscore_greater_than_70 += 1
								h_r_qscore_greater_than_50 += 1
								h_r_qscore_less_than_50 += 1
							elif(r_q_score >= .50):
								h_r_qscore_greater_than_50 += 1
								h_r_qscore_less_than_50 += 1
							elif(r_q_score >= .30):
								h_r_qscore_less_than_50 += 1
							elif(r_q_score < .3):
								h_r_qscore_less_than_30 += 1
						else:
							k_r_qscore.append(r_q_score)
							known_proteins_w_raptor += 1
							if(r_q_score >= .90):
								k_r_qscore_greater_than_90 += 1
								k_r_qscore_greater_than_70 += 1
								k_r_qscore_greater_than_50 += 1
								k_r_qscore_less_than_50 += 1
							elif(r_q_score >= .70):
								k_r_qscore_greater_than_70 += 1
								k_r_qscore_greater_than_50 += 1
								k_r_qscore_less_than_50 += 1
							elif(r_q_score >= .50):
								k_r_qscore_greater_than_50 += 1
								k_r_qscore_less_than_50 += 1
							elif(r_q_score >= .30):
								k_r_qscore_less_than_50 += 1
							elif(r_q_score < .3):
								k_r_qscore_less_than_30 += 1

								
	a_header = f'''
	Q-score distribution for AlphaFold Predictions
	'''

	a_header = a_header.replace("\t","")

	r_header = f'''
	Q-score distribution for RaptorX Predictions
	'''

	r_header = r_header.replace("\t","")

	# Plot text for hypothetical proteins
	h_atext = fr'''
	3D Homologs Found: {hypo_proteins_w_alpha}
	$\geq$0.90   {'{:4d}'.format(h_a_qscore_greater_than_90)} ({'{:5.2f}'.format(h_a_qscore_greater_than_90/len(h_a_qscore)*100)}%)
	$\geq$0.70   {'{:4d}'.format(h_a_qscore_greater_than_70)} ({'{:5.2f}'.format(h_a_qscore_greater_than_70/len(h_a_qscore)*100)}%)
	$\geq$0.50   {'{:4d}'.format(h_a_qscore_greater_than_50)} ({'{:5.2f}'.format(h_a_qscore_greater_than_50/len(h_a_qscore)*100)}%)
	$\geq$0.30   {'{:4d}'.format(h_a_qscore_less_than_50)} ({'{:5.2f}'.format((h_a_qscore_less_than_50/len(h_a_qscore)*100))}%)
	$<$0.30   {'{:4d}'.format(h_a_qscore_less_than_30)} ({'{:5.2f}'.format((h_a_qscore_less_than_30/len(h_a_qscore)*100))}%)
	'''

	h_atext = h_atext.replace("\t","")

	h_rtext = fr'''
	3D Homologs Found: {hypo_proteins_w_raptor}
	$\geq$0.90   {'{:4d}'.format(h_r_qscore_greater_than_90)} ({'{:5.2f}'.format(h_r_qscore_greater_than_90/len(h_r_qscore)*100)}%)
	$\geq$0.70   {'{:4d}'.format(h_r_qscore_greater_than_70)} ({'{:5.2f}'.format(h_r_qscore_greater_than_70/len(h_r_qscore)*100)}%)
	$\geq$0.50   {'{:4d}'.format(h_r_qscore_greater_than_50)} ({'{:5.2f}'.format(h_r_qscore_greater_than_50/len(h_r_qscore)*100)}%)
	$\geq$0.30   {'{:4d}'.format(h_r_qscore_less_than_50)} ({'{:5.2f}'.format((h_r_qscore_less_than_50/len(h_r_qscore)*100))}%)
	$<$0.30   {'{:4d}'.format(h_r_qscore_less_than_30)} ({'{:5.2f}'.format((h_r_qscore_less_than_30/len(h_r_qscore)*100))}%)
	'''

	h_rtext = h_rtext.replace("\t","")

	## Plot text for known proteins
	k_atext = fr'''
	3D Homologs Found: {known_proteins_w_alpha}
	$\geq$0.90   {'{:4d}'.format(k_a_qscore_greater_than_90)} ({'{:5.2f}'.format(k_a_qscore_greater_than_90/len(k_a_qscore)*100)}%)
	$\geq$0.70   {'{:4d}'.format(k_a_qscore_greater_than_70)} ({'{:5.2f}'.format(k_a_qscore_greater_than_70/len(k_a_qscore)*100)}%)
	$\geq$0.50   {'{:4d}'.format(k_a_qscore_greater_than_50)} ({'{:5.2f}'.format(k_a_qscore_greater_than_50/len(k_a_qscore)*100)}%)
	$\geq$0.30   {'{:4d}'.format(k_a_qscore_less_than_50)} ({'{:5.2f}'.format((k_a_qscore_less_than_50/len(k_a_qscore)*100))}%)
	$<$0.30   {'{:4d}'.format(k_a_qscore_less_than_30)} ({'{:5.2f}'.format((k_a_qscore_less_than_30/len(k_a_qscore)*100))}%)
	'''

	k_atext = k_atext.replace("\t","")

	k_rtext = fr'''
	3D Homologs Found: {known_proteins_w_raptor}
	$\geq$0.90   {'{:4d}'.format(k_r_qscore_greater_than_90)} ({'{:5.2f}'.format(k_r_qscore_greater_than_90/len(k_r_qscore)*100)}%)
	$\geq$0.70   {'{:4d}'.format(k_r_qscore_greater_than_70)} ({'{:5.2f}'.format(k_r_qscore_greater_than_70/len(k_r_qscore)*100)}%)
	$\geq$0.50   {'{:4d}'.format(k_r_qscore_greater_than_50)} ({'{:5.2f}'.format(k_r_qscore_greater_than_50/len(k_r_qscore)*100)}%)
	$\geq$0.30   {'{:4d}'.format(k_r_qscore_less_than_50)} ({'{:5.2f}'.format((k_r_qscore_less_than_50/len(k_r_qscore)*100))}%)
	$<$0.30   {'{:4d}'.format(k_r_qscore_less_than_30)} ({'{:5.2f}'.format((k_r_qscore_less_than_30/len(k_r_qscore)*100))}%)
	'''

	k_rtext = k_rtext.replace("\t","")

	# q score distribution with hypo-protein overlay
	bins = 50
	max = 1
	fig, (alp,rpx) = plt.subplots(2,1)

	blnk = []
	all = "#6F333B"
	hypo = "#1f77b4"
	known = "#ff7f0e"

	alp.hist([h_a_qscore,k_a_qscore],bins=50,rwidth = .95, color=[hypo,known], range = [0.1,1],stacked=True, label=[f"Hypothetical Proteins: {hypo_proteins}",f"Known Proteins: {known_proteins}"])
	alp.set_xlim(0.1,1)
	alp.set_ylim(0,100)
	alp.set_xticks([float(max/bins)*i for i in range(5,bins+1)])
	alp.text(.5,85,a_header,fontsize="x-large",fontfamily="monospace")
	alp.text(.5,50,h_atext,fontsize="x-large",fontfamily="monospace",color=hypo)
	alp.text(.7,50,k_atext,fontsize="x-large",fontfamily="monospace",color=known)
	alp.legend(loc=2)

	rpx.hist([h_r_qscore,k_r_qscore],bins=50, rwidth = .95, color=[hypo,known], range = [0.1,1],stacked=True, label=[f"Hypothetical Proteins: {hypo_proteins}",f"Known Proteins: {known_proteins}"])
	rpx.set_xlim(0.1,1)
	rpx.set_ylim(0,100)
	rpx.set_xticks([float(max/bins)*i for i in range(5,bins+1)])
	rpx.text(.5,85,r_header,fontsize="x-large",fontfamily="monospace")
	rpx.text(.5,50,h_rtext,fontsize="x-large",fontfamily="monospace",color=hypo)
	rpx.text(.7,50,k_rtext,fontsize="x-large",fontfamily="monospace",color=known)
	rpx.legend(loc=2)

	plt.show()

if(voro):
	###################################################################################################
	## Q-score distribution graph
	###################################################################################################

	a_vorocnn = []
	r_vorocnn = []

	all_proteins = 0

	a_vorocnn_greater_than_90 = 0
	a_vorocnn_greater_than_70 = 0
	a_vorocnn_greater_than_50 = 0
	a_vorocnn_less_than_50 = 0
	a_vorocnn_less_than_30 = 0

	r_vorocnn_greater_than_90 = 0
	r_vorocnn_greater_than_70 = 0
	r_vorocnn_greater_than_50 = 0
	r_vorocnn_less_than_50 = 0
	r_vorocnn_less_than_30 = 0

	## Known proteins
	k_a_vorocnn = []
	k_r_vorocnn = []

	known_proteins = 0
	known_proteins_w_alpha = 0
	known_proteins_w_raptor = 0

	k_a_vorocnn_greater_than_90 = 0
	k_a_vorocnn_greater_than_70 = 0
	k_a_vorocnn_greater_than_50 = 0
	k_a_vorocnn_less_than_50 = 0
	k_a_vorocnn_less_than_30 = 0

	k_r_vorocnn_greater_than_90 = 0
	k_r_vorocnn_greater_than_70 = 0
	k_r_vorocnn_greater_than_50 = 0
	k_r_vorocnn_less_than_50 = 0
	k_r_vorocnn_less_than_30 = 0

	## Hypothetical proteins
	h_a_vorocnn = []
	h_r_vorocnn = []

	hypo_proteins = 0
	hypo_proteins_w_alpha = 0
	hypo_proteins_w_raptor = 0

	h_a_vorocnn_greater_than_90 = 0
	h_a_vorocnn_greater_than_70 = 0
	h_a_vorocnn_greater_than_50 = 0
	h_a_vorocnn_less_than_50 = 0
	h_a_vorocnn_less_than_30 = 0

	h_r_vorocnn_greater_than_90 = 0
	h_r_vorocnn_greater_than_70 = 0
	h_r_vorocnn_greater_than_50 = 0
	h_r_vorocnn_less_than_50 = 0
	h_r_vorocnn_less_than_30 = 0

	with open (argv[2],"r") as FILE:
		for line in FILE:

			line = line.replace("\n",'')

			## Parse data for non-headers and divider lines
			if line != "" and line[0:5] != "Locus":
				data = line.split("\t")
				a_q_score = data[4]
				r_q_score = data[13]

				## Get data only for the best hits
				if(data[0] != ''):

					all_proteins += 1
					
					## Keep track of how many proteins are known vs unknown
					if(data[1] == "hypothetical protein"):
						hypo_proteins += 1
					elif("Uncharacterized" in data[1]):
						hypo_proteins += 1
					elif("HYPOTHETICAL" in data[1]):
						hypo_proteins += 1
					elif("UPF" in data[1]):
						hypo_proteins += 1
					else:
						known_proteins += 1
						

					## If there is a q score for AlphaFold, store it, and count what quality it falls under
					if(a_q_score != '---' and a_q_score != ''):

						a_q_score = float(data[4])
						a_vorocnn.append(a_q_score)
						
						if(a_q_score >= .90):
							a_vorocnn_greater_than_90 += 1
							a_vorocnn_greater_than_70 += 1
							a_vorocnn_greater_than_50 += 1
							a_vorocnn_less_than_50 += 1
						elif(a_q_score >= .70):
							a_vorocnn_greater_than_70 += 1
							a_vorocnn_greater_than_50 += 1
							a_vorocnn_less_than_50 += 1
						elif(a_q_score >= .50):
							a_vorocnn_greater_than_50 += 1
							a_vorocnn_less_than_50 += 1
						elif(a_q_score >= .30):
							a_vorocnn_less_than_50 += 1
						elif(a_q_score < .3):
							a_vorocnn_less_than_30 += 1
						
						if(data[1] == "hypothetical protein" or "Uncharacterized" in data[1] or "HYPOTHETICAL" in data[1] or "UPF" in data[1]):
							h_a_vorocnn.append(a_q_score)
							hypo_proteins_w_alpha += 1
							if(a_q_score >= .90):
								h_a_vorocnn_greater_than_90 += 1
								h_a_vorocnn_greater_than_70 += 1
								h_a_vorocnn_greater_than_50 += 1
								h_a_vorocnn_less_than_50 += 1
							elif(a_q_score >= .70):
								h_a_vorocnn_greater_than_70 += 1
								h_a_vorocnn_greater_than_50 += 1
								h_a_vorocnn_less_than_50 += 1
							elif(a_q_score >= .50):
								h_a_vorocnn_greater_than_50 += 1
								h_a_vorocnn_less_than_50 += 1
							elif(a_q_score >= .30):
								h_a_vorocnn_less_than_50 += 1
							elif(a_q_score < .3):
								h_a_vorocnn_less_than_30 += 1
						else:
							k_a_vorocnn.append(a_q_score)
							known_proteins_w_alpha += 1
							if(a_q_score >= .90):
								k_a_vorocnn_greater_than_90 += 1
								k_a_vorocnn_greater_than_70 += 1
								k_a_vorocnn_greater_than_50 += 1
								k_a_vorocnn_less_than_50 += 1
							elif(a_q_score >= .70):
								k_a_vorocnn_greater_than_70 += 1
								k_a_vorocnn_greater_than_50 += 1
								k_a_vorocnn_less_than_50 += 1
							elif(a_q_score >= .50):
								k_a_vorocnn_greater_than_50 += 1
								k_a_vorocnn_less_than_50 += 1
							elif(a_q_score >= .30):
								k_a_vorocnn_less_than_50 += 1
							elif(a_q_score < .3):
								k_a_vorocnn_less_than_30 += 1

					## If there is a q score for AlphaFold, store it, and count what quality it falls under
					if(r_q_score != '---' and r_q_score != ''):

						r_q_score = float(data[13])
						r_vorocnn.append(r_q_score)

						if(r_q_score >= .90):
							r_vorocnn_greater_than_90 += 1
							r_vorocnn_greater_than_70 += 1
							r_vorocnn_greater_than_50 += 1
							r_vorocnn_less_than_50 += 1
						elif(r_q_score >= .70):
							r_vorocnn_greater_than_70 += 1
							r_vorocnn_greater_than_50 += 1
							r_vorocnn_less_than_50 += 1
						elif(r_q_score >= .50):
							r_vorocnn_greater_than_50 += 1
							r_vorocnn_less_than_50 += 1
						elif(r_q_score >= .30):
							r_vorocnn_less_than_50 += 1
						elif(r_q_score < .3):
							r_vorocnn_less_than_30 += 1
						
						if(data[1] == "hypothetical protein" or "Uncharacterized" in data[1] or "HYPOTHETICAL" in data[1] or "UPF" in data[1]):
							h_r_vorocnn.append(r_q_score)
							hypo_proteins_w_raptor += 1
							if(r_q_score >= .90):
								h_r_vorocnn_greater_than_90 += 1
								h_r_vorocnn_greater_than_70 += 1
								h_r_vorocnn_greater_than_50 += 1
								h_r_vorocnn_less_than_50 += 1
							elif(r_q_score >= .70):
								h_r_vorocnn_greater_than_70 += 1
								h_r_vorocnn_greater_than_50 += 1
								h_r_vorocnn_less_than_50 += 1
							elif(r_q_score >= .50):
								h_r_vorocnn_greater_than_50 += 1
								h_r_vorocnn_less_than_50 += 1
							elif(r_q_score >= .30):
								h_r_vorocnn_less_than_50 += 1
							elif(r_q_score < .3):
								h_r_vorocnn_less_than_30 += 1
						else:
							k_r_vorocnn.append(r_q_score)
							known_proteins_w_raptor += 1
							if(r_q_score >= .90):
								k_r_vorocnn_greater_than_90 += 1
								k_r_vorocnn_greater_than_70 += 1
								k_r_vorocnn_greater_than_50 += 1
								k_r_vorocnn_less_than_50 += 1
							elif(r_q_score >= .70):
								k_r_vorocnn_greater_than_70 += 1
								k_r_vorocnn_greater_than_50 += 1
								k_r_vorocnn_less_than_50 += 1
							elif(r_q_score >= .50):
								k_r_vorocnn_greater_than_50 += 1
								k_r_vorocnn_less_than_50 += 1
							elif(r_q_score >= .30):
								k_r_vorocnn_less_than_50 += 1
							elif(r_q_score < .3):
								k_r_vorocnn_less_than_30 += 1

								
	a_header = f'''
	Q-score distribution for AlphaFold Predictions
	'''

	a_header = a_header.replace("\t","")

	r_header = f'''
	Q-score distribution for RaptorX Predictions
	'''

	r_header = r_header.replace("\t","")

	# Plot text for hypothetical proteins
	h_atext = fr'''
	3D Homologs Found: {hypo_proteins_w_alpha}
	$\geq$0.90   {'{:4d}'.format(h_a_vorocnn_greater_than_90)} ({'{:5.2f}'.format(h_a_vorocnn_greater_than_90/len(h_a_vorocnn)*100)}%)
	$\geq$0.70   {'{:4d}'.format(h_a_vorocnn_greater_than_70)} ({'{:5.2f}'.format(h_a_vorocnn_greater_than_70/len(h_a_vorocnn)*100)}%)
	$\geq$0.50   {'{:4d}'.format(h_a_vorocnn_greater_than_50)} ({'{:5.2f}'.format(h_a_vorocnn_greater_than_50/len(h_a_vorocnn)*100)}%)
	$\geq$0.30   {'{:4d}'.format(h_a_vorocnn_less_than_50)} ({'{:5.2f}'.format((h_a_vorocnn_less_than_50/len(h_a_vorocnn)*100))}%)
	$<$0.30   {'{:4d}'.format(h_a_vorocnn_less_than_30)} ({'{:5.2f}'.format((h_a_vorocnn_less_than_30/len(h_a_vorocnn)*100))}%)
	'''

	h_atext = h_atext.replace("\t","")

	h_rtext = fr'''
	3D Homologs Found: {hypo_proteins_w_raptor}
	$\geq$0.90   {'{:4d}'.format(h_r_vorocnn_greater_than_90)} ({'{:5.2f}'.format(h_r_vorocnn_greater_than_90/len(h_r_vorocnn)*100)}%)
	$\geq$0.70   {'{:4d}'.format(h_r_vorocnn_greater_than_70)} ({'{:5.2f}'.format(h_r_vorocnn_greater_than_70/len(h_r_vorocnn)*100)}%)
	$\geq$0.50   {'{:4d}'.format(h_r_vorocnn_greater_than_50)} ({'{:5.2f}'.format(h_r_vorocnn_greater_than_50/len(h_r_vorocnn)*100)}%)
	$\geq$0.30   {'{:4d}'.format(h_r_vorocnn_less_than_50)} ({'{:5.2f}'.format((h_r_vorocnn_less_than_50/len(h_r_vorocnn)*100))}%)
	$<$0.30   {'{:4d}'.format(h_r_vorocnn_less_than_30)} ({'{:5.2f}'.format((h_r_vorocnn_less_than_30/len(h_r_vorocnn)*100))}%)
	'''

	h_rtext = h_rtext.replace("\t","")

	## Plot text for known proteins
	k_atext = fr'''
	3D Homologs Found: {known_proteins_w_alpha}
	$\geq$0.90   {'{:4d}'.format(k_a_vorocnn_greater_than_90)} ({'{:5.2f}'.format(k_a_vorocnn_greater_than_90/len(k_a_vorocnn)*100)}%)
	$\geq$0.70   {'{:4d}'.format(k_a_vorocnn_greater_than_70)} ({'{:5.2f}'.format(k_a_vorocnn_greater_than_70/len(k_a_vorocnn)*100)}%)
	$\geq$0.50   {'{:4d}'.format(k_a_vorocnn_greater_than_50)} ({'{:5.2f}'.format(k_a_vorocnn_greater_than_50/len(k_a_vorocnn)*100)}%)
	$\geq$0.30   {'{:4d}'.format(k_a_vorocnn_less_than_50)} ({'{:5.2f}'.format((k_a_vorocnn_less_than_50/len(k_a_vorocnn)*100))}%)
	$<$0.30   {'{:4d}'.format(k_a_vorocnn_less_than_30)} ({'{:5.2f}'.format((k_a_vorocnn_less_than_30/len(k_a_vorocnn)*100))}%)
	'''

	k_atext = k_atext.replace("\t","")

	k_rtext = fr'''
	3D Homologs Found: {known_proteins_w_raptor}
	$\geq$0.90   {'{:4d}'.format(k_r_vorocnn_greater_than_90)} ({'{:5.2f}'.format(k_r_vorocnn_greater_than_90/len(k_r_vorocnn)*100)}%)
	$\geq$0.70   {'{:4d}'.format(k_r_vorocnn_greater_than_70)} ({'{:5.2f}'.format(k_r_vorocnn_greater_than_70/len(k_r_vorocnn)*100)}%)
	$\geq$0.50   {'{:4d}'.format(k_r_vorocnn_greater_than_50)} ({'{:5.2f}'.format(k_r_vorocnn_greater_than_50/len(k_r_vorocnn)*100)}%)
	$\geq$0.30   {'{:4d}'.format(k_r_vorocnn_less_than_50)} ({'{:5.2f}'.format((k_r_vorocnn_less_than_50/len(k_r_vorocnn)*100))}%)
	$<$0.30   {'{:4d}'.format(k_r_vorocnn_less_than_30)} ({'{:5.2f}'.format((k_r_vorocnn_less_than_30/len(k_r_vorocnn)*100))}%)
	'''

	k_rtext = k_rtext.replace("\t","")

	# q score distribution with hypo-protein overlay
	bins = 50
	max = 1
	fig, (alp,rpx) = plt.subplots(2,1)

	blnk = []
	all = "#6F333B"
	hypo = "#1f77b4"
	known = "#ff7f0e"

	alp.hist([h_a_vorocnn,k_a_vorocnn],bins=50,rwidth = .95, color=[hypo,known], range = [0.1,1],stacked=True, label=[f"Hypothetical Proteins: {hypo_proteins_w_alpha}",f"Known Proteins: {known_proteins_w_alpha}"])
	alp.set_xlim(0,1)
	alp.set_ylim(0,250)
	alp.set_xticks([float(max/bins)*i for i in range(5,bins+1)])
	alp.text(.12,200,a_header,fontsize="x-large",fontfamily="monospace")
	alp.text(.12,120,h_atext,fontsize="x-large",fontfamily="monospace",color=hypo)
	alp.text(.32,120,k_atext,fontsize="x-large",fontfamily="monospace",color=known)
	alp.legend(loc=2)

	rpx.hist([h_r_vorocnn,k_r_vorocnn],bins=50, rwidth = .95, color=[hypo,known], range = [0.1,1],stacked=True, label=[f"Hypothetical Proteins: {hypo_proteins_w_raptor}",f"Known Proteins: {known_proteins_w_raptor}"])
	rpx.set_xlim(0,1)
	rpx.set_ylim(0,250)
	rpx.set_xticks([float(max/bins)*i for i in range(5,bins+1)])
	rpx.text(.12,200,r_header,fontsize="x-large",fontfamily="monospace")
	rpx.text(.12,120,h_rtext,fontsize="x-large",fontfamily="monospace",color=hypo)
	rpx.text(.32,120,k_rtext,fontsize="x-large",fontfamily="monospace",color=known)
	rpx.legend(loc=2)

	plt.show()

if(tmvq):

	a_tm = []
	a_qs = []
	r_tm = []
	r_qs = []

	from numpy import polyfit

	def fit(x,curve):
		degree = len(curve)-1
		y = [0 for i in range(len(x))]
		for index,step in enumerate(x):
			store = 0
			for poly in range(degree,-1,-1):
				store += curve[degree-poly]*(step**poly)
			y[index] = store
		return y

	with open (argv[2],"r") as FILE:
		for line in FILE:

			line = line.replace("\n",'')

			## Parse data for non-headers and divider lines
			if line != "" and line[0:5] != "Locus":

				data = line.split("\t")

				a_qs_val = data[7]
				a_tm_val = data[8]

				if((a_tm_val != "---" and a_tm_val != "") and (a_qs_val != "---" and a_qs_val != "")):
					a_tm.append(float(a_tm_val))
					a_qs.append(float(a_qs_val))
	
	x = [i*.01 for i in range(0,101)]
	y = [i for i in x]

	curve = polyfit(a_qs,a_tm,deg=1)
	fitted = fit(x,curve)
	equation = r"y={0:.2g}x$^{1}$".format(curve[0],len(curve)-1)
	for i in range(1,len(curve)-1):
		equation += r"{0:+.2g}x$^{1}$".format(curve[i],len(curve)-i)
	equation += r"{0:+.2g}".format(curve[len(curve)-1])

	curve2 = polyfit(a_qs,a_tm,deg=2)
	fitted2 = fit(x,curve2)
	equation2 = r"y={0:.2g}x$^{1}$".format(curve2[0],len(curve2)-1)
	for i in range(1,len(curve2)-1):
		equation2 += r"{0:+.2g}x$^{1}$".format(curve2[i],len(curve2)-i)
	equation2 += r"{0:+.2g}".format(curve2[len(curve2)-1])

	fig, (a1) = plt.subplots(1,1)
	a1.scatter(a_qs,a_tm,s=10,color="#A53F56")
	a1.set_xlim(0,1)
	a1.set_ylim(0,1)
	a1.plot(x,y,color="blue",label="y=x")
	a1.plot(x,fitted,color="black",label = 'Fitted Curve (Linear): {0}'.format(equation))
	a1.plot(x,fitted2,color="deeppink",label = 'Fitted Curve (Quadratic): {0}'.format(equation2))
	a1.legend(loc=2)
	plt.show()

if(pdvvo):

	a_plddt = []
	a_voro = []

	with open (argv[2],"r") as FILE:
		for line in FILE:

			line = line.replace("\n",'')

			## Parse data for non-headers and divider lines
			if line != "" and line[0:5] != "Locus":

				data = line.split("\t")

				a_plddt_val = data[3]
				a_voro_val = data[4]

				if((a_plddt_val != "---" and a_plddt_val != "") and (a_voro_val != "---" and a_voro_val != "")):
					a_plddt.append(float(a_plddt_val))
					a_voro.append(float(a_voro_val))

	from numpy import polyfit

	def fit(x,curve):
		degree = len(curve)-1
		y = [0 for i in range(len(x))]
		for index,step in enumerate(x):
			store = 0
			for poly in range(degree,-1,-1):
				store += curve[degree-poly]*(step**poly)
			y[index] = store
		return y

	x = [i for i in range(0,101)]
	y = [.01*i for i in x]

	curve = polyfit(a_plddt,a_voro,deg=1)
	print(curve)
	fitted = fit(a_plddt,curve)
	equation = r"y={0:.2g}x$^{1}$".format(curve[0],len(curve)-1)
	for i in range(1,len(curve)-1):
		equation += r"{0:+.2g}x$^{1}$".format(curve[i],len(curve)-i)
	equation += r"{0:+.2g}".format(curve[len(curve)-1])

	fig, (a1) = plt.subplots(1,1)
	a1.scatter(a_plddt,a_voro,s=10,color="#A53F56")
	a1.plot(a_plddt,fitted,linestyle='-',label = 'Fitted Curve: {0}'.format(equation),color="black")
	a1.plot(x,y,color="blue",label="y=.01*x")
	a1.set_xlim(20,100)
	a1.set_ylim(0.2,.85)
	a1.legend(loc=2)
	plt.show()