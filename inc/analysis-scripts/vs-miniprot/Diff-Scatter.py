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
best_face_color = "#c2f08d"
best_edge_color = "#90dc39"
sum_face_color  = "#d2efb1"
sum_edge_color  = "#b9e38a"


species = sys.argv[2].lower()

if   (species == "human"       or species == "homo_sapiens"):
        
	# Lavender
	best_face_color = "#d19efb"
	best_edge_color = "#8327d0"
	sum_face_color  = "#e0c8f4"
	sum_edge_color  = "#b380dd"
	species         = "Human"
    
elif (species == "chicken"     or species == "gallus_gallus"):
	
	# Amber
	best_face_color = "#fcc06d"
	best_edge_color = "#ea9c12"
	sum_face_color  = "#f4e1a7"
	sum_edge_color  = "#e7c765"
	species         = "Chicken"


elif (species == "fruit_fly"   or species == "drosophila_melanogaster"):
	
	# Rose
	best_face_color = "#f8b8e6"
	best_edge_color = "#b61087"
	sum_face_color  = "#fccce2"
	sum_edge_color  = "#f391be"
	species         = "Fly"

elif (species == "stickleback" or species == "gasterosteus_aculeatus"):
	
	# Aqua
	best_face_color = "#2ae2ae"
	best_edge_color = "#3da98b"
	sum_face_color  = "#c4f3d7"
	sum_edge_color  = "#9ae6b9"
	species         = "Stickleback"




pct_regex = r"(\S+)\%$"


Data = pd.read_csv(sys.argv[1])


CoverageDiffs = []
PctIDDiffs    = []


for index,row in Data.iterrows():

	CoverageMatch = re.match(pct_regex,row['Peak-Coverage-Diff'])
	CoverageDiffs.append(float(CoverageMatch.group(1)))

	PctIDMatch = re.match(pct_regex,row['Pct-ID-Diff'])
	PctIDDiffs.append(float(PctIDMatch.group(1)))


plt.scatter(CoverageDiffs,PctIDDiffs,color=best_face_color,s=12)

plt.plot([-100,100],[   0,  0],color="black",linewidth=0.5)
plt.plot([   0,  0],[-100,100],color="black",linewidth=0.5)

plt.gca().set_xlim([-100.0,100.0])
plt.gca().set_ylim([-100.0,100.0])

plt.savefig(species+'.png')