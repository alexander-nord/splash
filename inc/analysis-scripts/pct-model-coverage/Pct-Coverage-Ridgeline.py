#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import re

import matplotlib as mpl
import matplotlib.pyplot as plt



if (len(sys.argv) != 3):
        sys.stderr.write("\n  USAGE: ./Pct-Coverage-Ridgeline.py [Coverage.csv] [species]\n\n")
        exit(1)


# Default
face_color = ""
edge_color = ""

species = sys.argv[2]
if   (species == "human"):
	# Lavender
	face_color = ""
	edge_color = ""
elif (species == "chicken"):
	# Amber
	face_color = ""
	edge_color = ""
elif (species == "fly"):
	# Rose
	face_color = ""
	edge_color = ""
elif (species == "cod"):
	# Aqua
	face_color = ""
	edge_color = ""


pct_regex = r"(\S+)\%"

Data = pd.read_csv(sys.argv[1])

BestESCoverage = []
SumESCoverage  = []




