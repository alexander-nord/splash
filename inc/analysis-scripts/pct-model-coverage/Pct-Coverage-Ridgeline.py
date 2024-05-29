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
	best_face_color = "#9df3db"
	best_edge_color = "#09e5a7"
	sum_face_color  = "#c4f3d7"
	sum_edge_color  = "#9ae6b9"
	species         = "Stickleback"



pct_regex = r"(\S+)\%"


Data = pd.read_csv(sys.argv[1])


BestESCoverage = []
SumESCoverage  = []


for index,row in Data.iterrows():

	PctMatch = re.match(pct_regex,row['Best-ES-Pct-Coverage'])
	BestESCoverage.append(float(PctMatch.group(1)))

	PctMatch = re.match(pct_regex,row['Sum-ES-Pct-Coverage'])
	SumESCoverage.append(float(PctMatch.group(1)))


fig = plt.figure()
gs  = (grid_spec.GridSpec(2,1))

AxObjs = []



#
# First is the best!
#
AxObjs.append(fig.add_subplot(gs[0:1,0:]))

BestData = pd.DataFrame(BestESCoverage,index=range(len(BestESCoverage)),columns=["Best-Model-Coverage"])
Plot     = (BestData.plot.kde(ax=AxObjs[-1],color=best_edge_color,lw=0.5))

X = Plot.get_children()[0]._x
Y = Plot.get_children()[0]._y

AxObjs[-1].fill_between(X,Y,color=best_face_color)
AxObjs[-1].set_xlim(0,100)
AxObjs[-1].set_ylim(0,0.05)

Rect = AxObjs[-1].patch
Rect.set_alpha(0)

AxObjs[-1].tick_params(left=False,bottom=False)
AxObjs[-1].set_xticklabels([])
AxObjs[-1].set_yticklabels([])
AxObjs[-1].set_ylabel('')

AxObjs[-1].spines["top"   ].set_visible(False)
AxObjs[-1].spines["bottom"].set_visible(False)
AxObjs[-1].spines["left"  ].set_visible(False)
AxObjs[-1].spines["right" ].set_visible(False)

AxObjs[-1].get_legend().remove()
# AxObjs[-1].text(-6,0.0005,"Best ECS Coverage",fontweight="bold",fontsize=14,ha="left")


#
# Second is the sum!
#
AxObjs.append(fig.add_subplot(gs[1:2,0:]))

SumData = pd.DataFrame(SumESCoverage,index=range(len(SumESCoverage)),columns=["Sum-Model-Coverage"])
Plot    = (SumData.plot.kde(ax=AxObjs[-1],color=sum_edge_color,lw=0.5))

X = Plot.get_children()[0]._x
Y = Plot.get_children()[0]._y

AxObjs[-1].fill_between(X,Y,color=sum_face_color)
AxObjs[-1].set_xlim(0,100)
AxObjs[-1].set_ylim(0,0.05)

Rect = AxObjs[-1].patch
Rect.set_alpha(0)

AxObjs[-1].tick_params(left=False)
AxObjs[-1].set_xlabel("Percent Model Coverage",fontsize=18,fontweight="bold")
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





