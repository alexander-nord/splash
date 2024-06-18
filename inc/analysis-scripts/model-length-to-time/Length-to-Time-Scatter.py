#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import re

import matplotlib as mpl
import matplotlib.pyplot as plt



if (len(sys.argv) != 3):
	sys.stderr.write("\n  USAGE: ./Length-to-Time-Scatter.py [Coverage.csv] [species]\n\n")
	exit(1)




# Default : Lime
bath_face_color   = "#90dc39"
splash_face_color = "#c2f08d"


species = sys.argv[2].lower()

if   (species == "human"       or species == "homo_sapiens"):
        
	# Lavender
	splash_face_color = "#d19efb"
	species = "Human"
    
elif (species == "chicken"     or species == "gallus_gallus"):
	
	# Amber
	splash_face_color = "#fcc06d"
	species = "Chicken"


elif (species == "fruit_fly"   or species == "drosophila_melanogaster"):
	
	# Rose
	splash_face_color = "#f8b8e6"
	species = "Fly"

elif (species == "stickleback" or species == "gasterosteus_aculeatus"):
	
	# Aqua
	splash_face_color = "#2ae2ae"
	species = "Stickleback"





time_regex = r"^(\d+):(\d+):(\d+\.\d+)$"


Data = pd.read_csv(sys.argv[1])

total_bath_seconds = 0.0;
total_splash_seconds = 0.0;


ModelLengths  = []
BathSeconds   = []
SplashSeconds = []


for index,row in Data.iterrows():


	bath_time = row['BATH-Runtime']
	BathMatch = re.match(time_regex,bath_time)

	bath_seconds = 3600 * int(BathMatch.group(1)) + 60 * int(BathMatch.group(2)) + float(BathMatch.group(3))
	total_bath_seconds += bath_seconds

	BathSeconds.append(bath_seconds)


	#if (bath_seconds > 1000):
	#print(f"{row['Gene']}:{row['Sequence-ID']}")


	splash_time = row['Time-in-Splash']
	SplashMatch = re.match(time_regex,splash_time)

	splash_seconds = 3600 * int(SplashMatch.group(1)) + 60 * int(SplashMatch.group(2)) + float(SplashMatch.group(3))
	total_splash_seconds += splash_seconds

	SplashSeconds.append(splash_seconds)


	ModelLengths.append(row['Model-Length'])



total_bath_seconds   = int(10 * total_bath_seconds  ) / 10.0
total_splash_seconds = int(10 * total_splash_seconds) / 10.0

print(f" BATH  : {total_bath_seconds  } seconds")
print(f"splash : {total_splash_seconds} seconds")


plt.scatter(ModelLengths,BathSeconds,s=12,color=bath_face_color)
plt.scatter(ModelLengths,SplashSeconds,s=12,color=splash_face_color)
plt.xscale('log')
#plt.ylim((0,20))

plt.savefig(species+'.png')




