#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import re

import matplotlib as mpl
import matplotlib.pyplot as plt



if (len(sys.argv) != 2):
	sys.stderr.write("\n  USAGE: ./Length-to-Time-Scatter.py [Coverage.csv]\n\n")
	exit(1)



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


	if (bath_seconds > 100):
		print(f"{row['Gene']}:{row['Sequence-ID']}")


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


plt.scatter(ModelLengths,BathSeconds)
plt.scatter(ModelLengths,SplashSeconds)
plt.xscale('log')
plt.ylim((0,20))
plt.show()




