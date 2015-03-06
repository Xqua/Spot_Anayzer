#!/usr/bin/env python 

import SpotAnaliser
import sys
import os
import matplotlib.pyplot as plt
import prettyplotlib as ppl
import numpy as np

if len(sys.argv) < 3:
	print """USAGE: train.py PATH_TO_IMAGE_FILES TXT_IMAGE_FILE_LIST
	to run the example:
	python train.py example/images list.txt
	the file list.txt needs to have for each line: image_name,number of filament"""

path = os.path.join(os.getcwd(),sys.argv[1])
filelist = sys.argv[2]


os.chdir(path)

lines = open(filelist).readlines()
bacterial_files = [i.strip().split(',')[0] for i in lines if (i and 'tif' in i.lower())]
known_filament_number = [float(i.strip().split(',')[1]) for i in lines if (i and 'tif' in i.lower())]

def loss_function(MreB, MreB_True):
	E = 0
	for i in range(len(MreB)):
		E += ((MreB[i] - MreB_True[i])**2)/len(MreB)
	return E

# self.timethreshold = 4 *self.timeinterval*(1/self.deltax) #Defined how far appart can peak be (in second)
# self.Time_Error = 6 / self.timeinterval #Allowed variance in second for a peak to be called true when found in vertical axis and checked VS time  
# self.y_threshold = 0.05				#Minimum signal intensity to be called a peak (the signal is normalized and goes from 0-1 should be set to mean(noise) + 2std(noise) or learned )
# self.dy_threshold = 0.001	
timeinterval = 2
deltax = 0.1
# timethreshold_totest = np.arange(1,10,1) * timeinterval * (1/deltax)
# timererror_totest = np.arange(1,15) / timeinterval
# y_threshold_totest = np.arange(0,0.4,0.01)
# dy_threshold_totest = np.arange(0,0.01,0.0001)
timethreshold_totest = np.arange(1,5,0.5) * timeinterval * (1/deltax)
timererror_totest = np.arange(1,7) / timeinterval
y_threshold_totest = np.arange(0,0.4,0.01)
dy_threshold_totest = np.arange(0,0.001,0.0001)

Error = np.zeros((len(timethreshold_totest),len(timererror_totest), len(y_threshold_totest), len(dy_threshold_totest)))

for i in range(len(timethreshold_totest)):
	for j in range(len(timererror_totest)):
		for k in range(len(y_threshold_totest)):
			for l in range(len(dy_threshold_totest)):

				MreB_mu = []
				MreB_std = []
				lenght = []
				speed = []

				for filename in bacterial_files:
					A = SpotAnaliser.Spot_Analysis(filename)
					A.timethreshold = timethreshold_totest[i]
					A.Time_Error = timererror_totest[j]
					A.y_threshold = y_threshold_totest[k]
					A.dy_threshold = dy_threshold_totest[l]

					A.Main()
					MreB_mu.append(A.meanMREB)
					# MreB_std.append(A.stdMREB)
					# speed.append(A.speed)
					# lenght.append(A.bact_lenght_nm)
				E = loss_function(MreB_mu, known_filament_number)
				Error[i][j][k][l] = E

# print Error
# print np.where(Error == Error.min())
tmp = np.where(Error == Error.min())
i,j,k,l = tmp[0][0], tmp[1][0], tmp[2][0], tmp[3][0]
print "timethreshold = ", (timethreshold_totest[i] / timeinterval) / (1/deltax)
print "Time_Error = ", timererror_totest[j] * timeinterval
print "y_threshold = ", y_threshold_totest[k]
print "dy_threshold = ", dy_threshold_totest[l]