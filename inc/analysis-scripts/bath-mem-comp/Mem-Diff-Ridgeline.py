#!/usr/bin/env python3

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as grid_spec

import pandas as pd
import numpy as np

import sys
import re



if (len(sys.argv) != 3):
        sys.stderr.write("\n  USAGE: ./Pct-Coverage-Ridgeline.py [Coverage.csv] [species]\n\n")
        exit(1)


y_lim = 0.006


# Default : Lime
face_color = "#c2f08d"
edge_color = "#90dc39"


species = sys.argv[2].lower()

if   (species == "human"       or species == "homo_sapiens"):
	
	# Lavender
	face_color = "#d19efb"
	edge_color = "#8327d0"
	species    = "Human"
	

elif (species == "chicken"     or species == "gallus_gallus"):
	
	# Amber
	face_color = "#fcc06d"
	edge_color = "#ea9c12"
	species    = "Chicken"


elif (species == "fruit_fly"   or species == "drosophila_melanogaster"):
	
	# Rose
	face_color = "#f8b8e6"
	edge_color = "#b61087"
	species    = "Fly"


elif (species == "stickleback" or species == "gasterosteus_aculeatus"):
	
	# Aqua
	face_color = "#2ae2ae"
	edge_color = "#3da98b"
	species    = "Stickleback"



Data = pd.read_csv(sys.argv[1])


MemDiffs = []
for index,row in Data.iterrows():
	MemDiffs.append(float(row['diff']))



fig = plt.figure(figsize=(9,5))
gs  = (grid_spec.GridSpec(2,1))

AxObjs = []



#
# Second is the sum!
#
AxObjs.append(fig.add_subplot(gs[0:1,0:]))

DiffData = pd.DataFrame(MemDiffs,index=range(len(MemDiffs)),columns=["diff"])
Plot     = (DiffData.plot.kde(ax=AxObjs[-1],color=edge_color,lw=0.5))

X = Plot.get_children()[0]._x
Y = Plot.get_children()[0]._y

AxObjs[-1].fill_between(X,Y,color=face_color)

AxObjs[-1].set_xlim(-1000,1000)
AxObjs[-1].set_ylim(0,y_lim)


Rect = AxObjs[-1].patch
Rect.set_alpha(0)

AxObjs[-1].tick_params(left=False)
AxObjs[-1].set_xlabel("Memory Usage Difference (MB)",fontsize=18,fontweight="bold")
AxObjs[-1].set_yticklabels([])
AxObjs[-1].set_ylabel('')

AxObjs[-1].spines["top"   ].set_visible(False)
AxObjs[-1].spines["bottom"].set_visible(False)
AxObjs[-1].spines["left"  ].set_visible(False)
AxObjs[-1].spines["right" ].set_visible(False)

AxObjs[-1].get_legend().remove()
# AxObjs[-1].text(-6,0.0005,"Sum ECS Coverage",fontweight="bold",fontsize=14,ha="left")


#
# Print it!
#
gs.update(hspace=-0.5)
plt.savefig(species+'.png')





