#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import re

import matplotlib as mpl
import matplotlib.pyplot as plt



if (len(sys.argv) != 3):
	sys.stderr.write("\n  USAGE: ./Length-to-Coverage-Scatter.py [species-differences.csv] [species]\n\n")
	exit(1)




# Default : Lime
species_color = "#c2f08d"
species       = "None"

species = sys.argv[2].lower()

if   (species == "human"       or species == "homo_sapiens"):
        
	# Lavender
	species_color = "#d19efb"
	species         = "Human"
    

elif (species == "chicken"     or species == "gallus_gallus"):
	
	# Amber
	species_color = "#fcc06d"
	species       = "Chicken"


elif (species == "fruit_fly"   or species == "drosophila_melanogaster"):
	
	# Rose
	species_color = "#f8b8e6"
	species       = "Fly"


elif (species == "stickleback" or species == "gasterosteus_aculeatus"):
	
	# Aqua
	species_color = "#2ae2ae"
	species       = "Stickleback"




pct_regex = r"(\S+)\%$"


Data = pd.read_csv(sys.argv[1])


CoverageDiffs = []
PctIDDiffs    = []


pos_coverage_pos_pct = 0;
pos_coverage_nil_pct = 0;
pos_coverage_neg_pct = 0;

nil_coverage_pos_pct = 0;
nil_coverage_nil_pct = 0;
nil_coverage_neg_pct = 0;

neg_coverage_pos_pct = 0;
neg_coverage_nil_pct = 0;
neg_coverage_neg_pct = 0;


for index,row in Data.iterrows():


	CoverageMatch = re.match(pct_regex,row['Peak-Coverage-Diff'])
	coverage_diff = float(CoverageMatch.group(1))

	PctIDMatch  = re.match(pct_regex,row['Pct-ID-Diff'])
	pct_id_diff = float(PctIDMatch.group(1))

	CoverageDiffs.append(coverage_diff)
	PctIDDiffs.append(pct_id_diff)


	if (coverage_diff > 0.0):
		if   (pct_id_diff > 0.0):
			pos_coverage_pos_pct += 1
		elif (pct_id_diff < 0.0):
			pos_coverage_neg_pct += 1
		else:
			pos_coverage_nil_pct += 1
	elif (coverage_diff < 0.0):
		if   (pct_id_diff > 0.0):
			neg_coverage_pos_pct += 1
		elif (pct_id_diff < 0.0):
			neg_coverage_neg_pct += 1
		else:
			neg_coverage_nil_pct += 1
	else:
		if   (pct_id_diff > 0.0):
			nil_coverage_pos_pct += 1
		elif (pct_id_diff < 0.0):
			nil_coverage_neg_pct += 1
		else:
			nil_coverage_nil_pct += 1


plt.scatter(CoverageDiffs,PctIDDiffs,color=species_color,s=12)

plt.plot([-100,100],[   0,  0],color="black",linewidth=0.5)
plt.plot([   0,  0],[-100,100],color="black",linewidth=0.5)

plt.gca().set_xlim([-100.0,100.0])
plt.gca().set_ylim([-100.0,100.0])

plt.savefig(species+'.png')


print("\n")
print( " Cov / Pct")
print("-");
print(f"  +  /  +  : {pos_coverage_pos_pct}")
print(f"  +  /  0  : {pos_coverage_nil_pct}")
print(f"  +  /  -  : {pos_coverage_neg_pct}")
print("-");
print(f"  0  /  +  : {nil_coverage_pos_pct}")
print(f"  0  /  0  : {nil_coverage_nil_pct}")
print(f"  0  /  -  : {nil_coverage_neg_pct}")
print("-");
print(f"  -  /  +  : {neg_coverage_pos_pct}")
print(f"  -  /  0  : {neg_coverage_nil_pct}")
print(f"  -  /  -  : {neg_coverage_neg_pct}")
print("\n")

